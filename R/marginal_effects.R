# Marginal cML-NMR standardization -------------------------------------------

#' Marginal effects from a component ML-NMR fit
#'
#' Standardizes treatment-specific posterior outcomes over an explicit target
#' distribution, then forms treatment contrasts within each posterior draw.
#' This is posterior standardization, not evaluation at a single covariate
#' profile. Treatment-specific outcomes are averaged first, and contrasts are
#' calculated second.
#'
#' @section Measures:
#' Available `measure` values are:
#'
#' * binomial: `"odds_ratio"`, `"risk_ratio"`, and `"risk_difference"`;
#' * Gaussian: `"mean_difference"`;
#' * Poisson: `"rate_ratio"` and `"rate_difference"`, for unit exposure;
#' * survival: `"survival_difference"`, `"survival_ratio"`,
#'   `"risk_difference"`, `"risk_ratio"`, `"rmst_difference"`,
#'   `"rmst_ratio"`, and `"time_specific_hazard_ratio"`.
#'
#' A scalar marginal hazard ratio is not defined. The time-specific marginal
#' hazard ratio generally changes with time even when the fitted conditional
#' hazards are proportional.
#'
#' Survival and risk ratios, restricted mean survival time ratios, and
#' time-specific marginal hazard ratios are accumulated and contrasted on the
#' log scale with log-sum-exp calculations before optional back-transformation.
#' Set `backtransf = FALSE` to retain the log-ratio reporting scale.
#' RMST measures require a piecewise-exponential donor baseline and are
#' evaluated analytically over its intervals. M-spline fits support marginal
#' survival, risk, and time-specific marginal hazard measures, but not RMST.
#'
#' @section Target distribution:
#' A data-frame target supplies empirical integration rows. `weights` are
#' normalized, and zero-weight rows are removed. Bernoulli modifiers must be
#' represented by actual 0/1 rows; a fractional prevalence is not an empirical
#' pseudo-patient. Continuous modifiers are accepted at their observed values,
#' subject to the support of the fitted margin.
#'
#' A list target defines a Gaussian-copula distribution through named `means`,
#' `sds`, and `margins`. Optional `cor` is a latent-scale correlation matrix and
#' `n_int` is the deterministic Sobol' integration size. If `cor` is omitted,
#' the correlation retained by [cmlnmr()] is used. Each continuous modifier may
#' use a `normal`, `gamma`, `lognormal`, or `beta` margin; fitted Bernoulli
#' modifiers must remain Bernoulli. Independent Bernoulli modifiers are
#' enumerated exactly, and continuous modifiers use deterministic Sobol' nodes.
#' Summary nodes must reproduce requested Bernoulli strata and non-normal means
#' and SDs within the implementation's moment-fidelity tolerances: Bernoulli
#' prevalence is checked to `max(1e-8, 0.1 * min(p, 1 - p))`, non-normal means
#' to 0.1 target SD, and non-normal SDs to 15 percent relative error. A failed
#' gate stops with a request for a larger `n_int` or an empirical target. Mixed
#' or non-normal multivariable targets also require resolved pairwise
#' correlations to agree within 0.05 of a high-resolution deterministic
#' reference grid. The resolved nodes and weights are retained on the result.
#' Target-distribution uncertainty is not propagated.
#'
#' @section Baseline transport:
#' Binomial measures, Poisson rate differences, and every survival measure need
#' `baseline_study`. Its fitted intercept is transported to the target; survival
#' measures also transport its fitted baseline hazard. The target covariate
#' distribution is always supplied separately. This is an explicit transport
#' assumption, not something identified by relative treatment effects.
#' Gaussian mean differences and Poisson rate ratios do not need a baseline
#' study because the common intercept cancels.
#'
#' Survival predictions start at model time zero. They are not landmark
#' predictions and are not conditioned on delayed entry or a supplied row-level
#' entry time. Prediction times must be positive and within the observed
#' follow-up support retained for the selected donor study, which can be shorter
#' than the global survival support.
#'
#' @section Identification and random effects:
#' Nonlinear marginal contrasts are screened at every positive-weight target
#' node. A passing result is labelled `"first-order screen"`, never `"exact"`,
#' because aggregate-data interactions and survival baselines can still be
#' weakly identified. Gaussian identity-link mean differences can be labelled
#' `"exact"` when their target-mean contrast is supported by IPD. Nonlinear
#' results that fail this conservative treatment-surface check are returned as
#' `NA` with `basis = "first-order screen failed"`. This screen covers the
#' treatment `beta`/`Gamma` surface, not full identification of prognostic
#' effects, donor intercepts, or donor baseline-hazard parameters.
#'
#' The target checks are support and moment checks, not a formal overlap or
#' positivity diagnostic. No overlap statistic is estimated, and extrapolation
#' risk remains a substantive limitation for the analyst to assess.
#'
#' Only `random_effect = "population"` is available. It sets study-arm
#' deviations to zero. The function does not make predictive draws over a new
#' study's heterogeneity distribution and does not return marginal component
#' effects or marginal rankings, since nonlinear standardization does not
#' preserve component additivity. Nonlinear marginal effects remain
#' treatment-level contrasts.
#'
#' @param object A fitted [cmlnmr()] object.
#' @param target Either a data frame containing exactly the fitted effect
#'   modifiers, or a list with named `means`, `sds`, `margins`, and optional
#'   `cor` and `n_int` entries.
#' @param weights Optional nonnegative target-row weights. This is only valid
#'   when `target` is a data frame.
#' @param reference Reference treatment.
#' @param all_contrasts Return every ordered treatment contrast if `TRUE`.
#' @param measure Marginal contrast measure. Available measures depend on the
#'   outcome family.
#' @param baseline_study Study whose fitted intercept, and for survival whose
#'   fitted baseline hazard, is transported to the target population.
#' @param times Positive prediction times for survival measures, restricted to
#'   the observed follow-up support of `baseline_study`. Predictions begin at
#'   model time zero rather than a landmark or delayed-entry time. For RMST
#'   these are the integration horizons. RMST measures require a
#'   piecewise-exponential donor baseline and are analytic over its intervals.
#' @param random_effect Random-effect prediction policy. Only the population
#'   mean, which sets study-arm deviations to zero, is currently available.
#' @param backtransf Report ratio measures on their natural scale if `TRUE`.
#' @param level Credible interval level.
#' @param ... Unused.
#'
#' @return A `cpaic_effects` data frame. `estimate`, `lower`, and `upper` are on
#'   the reporting scale. `estimate_contrast`, `se_contrast`, and
#'   `contrast_scale` describe the draw-level contrast used for inference: a log
#'   contrast for ratios and a natural-scale difference for difference measures.
#'   Attributes record the measure, target nodes and weights, donor baseline,
#'   random-effect policy, and sampler diagnostic status.
#' @export
marginal_effects <- function(
    object, target, weights = NULL, reference = NULL,
    all_contrasts = FALSE, measure = NULL, baseline_study = NULL,
    times = NULL, random_effect = "population", backtransf = TRUE,
    level = 0.95, ...) {
  UseMethod("marginal_effects")
}

#' @export
marginal_effects.default <- function(
    object, target, weights = NULL, reference = NULL,
    all_contrasts = FALSE, measure = NULL, baseline_study = NULL,
    times = NULL, random_effect = "population", backtransf = TRUE,
    level = 0.95, ...) {
  stop("`marginal_effects()` requires a cmlnmr fit.", call. = FALSE)
}

.cpaic_marginal_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.cpaic_marginal_calibrate_normal_nodes <- function(
    X, weights, means, sds, margins, cor) {
  normal <- which(margins == "normal")
  if (!length(normal)) return(X)
  if (nrow(X) < 2L) {
    stop("A summary target with a nondegenerate normal margin needs at least ",
         "two resolved integration nodes. Increase `n_int`.", call. = FALSE)
  }

  for (j in normal) {
    center <- sum(weights * X[, j])
    spread <- sqrt(sum(weights * (X[, j] - center)^2))
    if (!is.finite(spread) || spread <= 0) {
      stop("The target integration grid cannot represent the requested ",
           "normal margin. Increase `n_int`.", call. = FALSE)
    }
    X[, j] <- means[j] + sds[j] * (X[, j] - center) / spread
  }

  if (length(normal) == ncol(X) && length(normal) > 1L) {
    centered <- sweep(X, 2L, means, "-")
    standardized <- sweep(centered, 2L, sds, "/")
    weighted <- standardized * sqrt(weights)
    sample_cor <- crossprod(weighted)
    sample_chol <- tryCatch(chol(sample_cor), error = function(e) NULL)
    if (is.null(sample_chol)) {
      stop("The target integration grid cannot represent the requested ",
           "normal correlation. Increase `n_int`.", call. = FALSE)
    }
    target_cor <- cor %||% diag(length(normal))
    standardized <- standardized %*% solve(sample_chol) %*%
      chol(target_cor)
    X <- sweep(sweep(standardized, 2L, sds, "*"), 2L, means, "+")
  }
  X
}

.cpaic_marginal_weighted_cor <- function(X, weights) {
  means <- as.numeric(crossprod(weights, X))
  centered <- sweep(X, 2L, means, "-")
  covariance <- crossprod(centered * sqrt(weights))
  scales <- sqrt(diag(covariance))
  out <- covariance / outer(scales, scales)
  out[!is.finite(out)] <- NA_real_
  diag(out)[scales > 0] <- 1
  dimnames(out) <- list(colnames(X), colnames(X))
  out
}

.cpaic_marginal_joint_fidelity <- function(
    X, weights, means, sds, margins, cor, n_int) {
  Q <- ncol(X)
  if (Q < 2L || all(margins == "normal")) {
    return(list(resolved = NULL, reference = NULL))
  }
  n_reference <- as.integer(
    max(8192, min(65536, as.double(n_int) * 8))
  )
  target_cor <- cor %||% diag(Q)
  reference <- .cpaic_integration_points(
    means, sds, n_reference, cor = target_cor, margins = margins
  )
  reference_weights <- rep(1 / n_reference, n_reference)
  reference <- .cpaic_marginal_calibrate_normal_nodes(
    reference, reference_weights, means, sds, margins, target_cor
  )
  colnames(reference) <- names(means)
  resolved_cor <- .cpaic_marginal_weighted_cor(X, weights)
  reference_cor <- .cpaic_marginal_weighted_cor(
    reference, reference_weights
  )
  compare <- upper.tri(resolved_cor) & is.finite(resolved_cor) &
    is.finite(reference_cor)
  if (any(compare) &&
      any(abs(resolved_cor[compare] - reference_cor[compare]) > 0.05)) {
    stop("The target integration grid does not reproduce the requested joint ",
         "dependence closely enough. Increase `n_int` or supply a weighted ",
         "empirical target.", call. = FALSE)
  }
  list(resolved = resolved_cor, reference = reference_cor)
}

.cpaic_marginal_summary_nodes <- function(
    object, means, sds, margins, cor, n_int) {
  Q <- length(means)
  if (object$family == "gaussian") {
    return(list(
      X = matrix(means, nrow = 1L,
                 dimnames = list(NULL, names(means))),
      weights = 1,
      integration = "analytic target mean"
    ))
  }

  independent <- is.null(cor) ||
    max(abs(cor - diag(Q))) <= sqrt(.Machine$double.eps)
  bernoulli <- which(margins == "bernoulli")
  other <- setdiff(seq_len(Q), bernoulli)

  if (length(bernoulli) && independent) {
    n_binary <- 2^length(bernoulli)
    n_other <- if (length(other)) n_int else 1L
    if (!is.finite(n_binary) || n_binary * n_other > 131072) {
      stop("Exact Bernoulli stratification would create too many target ",
           "nodes. Supply an empirical target or reduce its dimension.",
           call. = FALSE)
    }
    binary_grid <- as.matrix(expand.grid(
      rep(list(c(0, 1)), length(bernoulli)), KEEP.OUT.ATTRS = FALSE
    ))
    binary_weights <- apply(binary_grid, 1L, function(z) {
      prod(ifelse(z == 1, means[bernoulli], 1 - means[bernoulli]))
    })
    if (length(other)) {
      X_other <- .cpaic_integration_points(
        means[other], sds[other], n_int,
        cor = if (length(other) > 1L) diag(length(other)) else NULL,
        margins = margins[other]
      )
    } else {
      X_other <- matrix(numeric(), nrow = 1L, ncol = 0L)
    }
    X <- matrix(NA_real_, nrow = nrow(binary_grid) * nrow(X_other), ncol = Q)
    w <- numeric(nrow(X))
    cursor <- 1L
    for (b in seq_len(nrow(binary_grid))) {
      rows <- cursor:(cursor + nrow(X_other) - 1L)
      if (length(other)) X[rows, other] <- X_other
      X[rows, bernoulli] <- matrix(
        binary_grid[b, ], nrow = length(rows), ncol = length(bernoulli),
        byrow = TRUE
      )
      w[rows] <- binary_weights[b] / nrow(X_other)
      cursor <- max(rows) + 1L
    }
    integration <- "exact independent Bernoulli strata with Sobol nodes"
  } else {
    X <- .cpaic_integration_points(
      means, sds, n_int, cor = cor, margins = margins
    )
    w <- rep(1 / nrow(X), nrow(X))
    integration <- "Sobol Gaussian-copula nodes"
  }

  X <- .cpaic_marginal_calibrate_normal_nodes(
    X, w, means, sds, margins, cor
  )
  colnames(X) <- names(means)
  joint_fidelity <- .cpaic_marginal_joint_fidelity(
    X, w, means, sds, margins, cor, n_int
  )

  realized_mean <- as.numeric(crossprod(w, X))
  realized_sd <- vapply(seq_len(Q), function(j) {
    sqrt(sum(w * (X[, j] - realized_mean[j])^2))
  }, numeric(1))
  names(realized_mean) <- names(realized_sd) <- names(means)

  for (j in seq_len(Q)) {
    if (margins[j] == "bernoulli") {
      p <- means[j]
      support <- sort(unique(X[w > 0, j]))
      if (p > 0 && p < 1 && !identical(support, c(0, 1))) {
        stop("The target integration grid omits a positive-probability ",
             "Bernoulli stratum for `", names(means)[j], "`. Increase ",
             "`n_int` or supply weighted empirical 0/1 rows.", call. = FALSE)
      }
      tolerance <- max(1e-8, 0.1 * min(p, 1 - p))
      if (abs(realized_mean[j] - p) > tolerance) {
        stop("The target integration grid does not reproduce the requested ",
             "Bernoulli prevalence for `", names(means)[j], "`. Increase ",
             "`n_int` or supply weighted empirical 0/1 rows.", call. = FALSE)
      }
    } else if (margins[j] != "normal") {
      mean_error <- abs(realized_mean[j] - means[j]) / sds[j]
      sd_error <- abs(realized_sd[j] - sds[j]) / sds[j]
      if (!is.finite(mean_error) || !is.finite(sd_error) ||
          mean_error > 0.1 || sd_error > 0.15) {
        stop("The target integration grid does not reproduce the requested ",
             "moments for `", names(means)[j], "`. Increase `n_int` or ",
             "supply an empirical target.", call. = FALSE)
      }
    }
  }

  list(X = X, weights = w, integration = integration,
       realized_means = realized_mean, realized_sds = realized_sd,
       realized_cor = joint_fidelity$resolved,
       reference_cor = joint_fidelity$reference)
}

.cpaic_marginal_target <- function(object, target, weights) {
  ems <- object$effect_modifiers
  margins <- object$margins[ems]
  if (is.null(target)) {
    stop("`target` must supply an explicit target covariate distribution.",
         call. = FALSE)
  }

  if (is.list(target) && !is.data.frame(target)) {
    if (!is.null(weights)) {
      stop("`weights` can only be supplied with an empirical data-frame target.",
           call. = FALSE)
    }
    allowed <- c("means", "sds", "margins", "cor", "n_int")
    unknown <- setdiff(names(target), allowed)
    if (length(unknown)) {
      stop("Summary `target` has unknown field(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    if (is.null(target$means) || is.null(target$sds)) {
      stop("Summary `target` must contain named `means` and `sds`.",
           call. = FALSE)
    }
    means <- target$means
    sds <- target$sds
    tmargins <- target$margins %||% margins
    for (x in list(means = means, sds = sds, margins = tmargins)) {
      if (is.null(names(x)) || !setequal(names(x), ems)) {
        stop("Summary target fields must be named for exactly the fitted effect ",
             "modifiers: ", paste(ems, collapse = ", "), ".", call. = FALSE)
      }
    }
    means <- as.numeric(means[ems]); names(means) <- ems
    sds <- as.numeric(sds[ems]); names(sds) <- ems
    tmargins <- as.character(tmargins[ems]); names(tmargins) <- ems
    if (any(!is.finite(means)) || any(!is.finite(sds))) {
      stop("Summary target means and SDs must be finite.", call. = FALSE)
    }
    supported_margins <- c("normal", "bernoulli", "gamma", "lognormal",
                           "beta")
    bad_margin <- setdiff(tmargins, supported_margins)
    if (length(bad_margin)) {
      stop("Unsupported summary target margin(s): ",
           paste(unique(bad_margin), collapse = ", "), ".", call. = FALSE)
    }
    if (any((tmargins == "bernoulli") != (margins == "bernoulli"))) {
      stop("Summary target margins must preserve which fitted effect ",
           "modifiers are Bernoulli.", call. = FALSE)
    }
    nonbern <- tmargins != "bernoulli"
    if (any(sds[nonbern] <= 0) || any(sds[!nonbern] < 0)) {
      stop("Summary target SDs must be positive for non-Bernoulli margins.",
           call. = FALSE)
    }
    if (any(tmargins == "bernoulli" & (means < 0 | means > 1))) {
      stop("Bernoulli target means must lie in [0, 1].", call. = FALSE)
    }
    if (any(tmargins %in% c("gamma", "lognormal") & means <= 0)) {
      stop("Gamma and lognormal target means must be positive.", call. = FALSE)
    }
    if (any(tmargins == "beta" & (means <= 0 | means >= 1))) {
      stop("Beta target means must lie in (0, 1).", call. = FALSE)
    }
    beta_idx <- tmargins == "beta"
    if (any(beta_idx & sds^2 >= means * (1 - means))) {
      stop("Beta target variances must be smaller than mean * (1 - mean).",
           call. = FALSE)
    }
    n_int <- target$n_int %||% 64L
    n_int_integer <- suppressWarnings(as.integer(n_int))
    if (!is.numeric(n_int) || length(n_int) != 1L || !is.finite(n_int) ||
        n_int < 1 || is.na(n_int_integer) || n_int != n_int_integer) {
      stop("Summary target `n_int` must be a positive integer.", call. = FALSE)
    }
    cor <- target$cor %||% object$cor
    if (!is.null(cor)) {
      cor <- as.matrix(cor)
      rn <- rownames(cor)
      cn <- colnames(cor)
      if (!is.null(rn) || !is.null(cn)) {
        if (is.null(rn) || is.null(cn) || !setequal(rn, ems) ||
            !setequal(cn, ems)) {
          stop("Summary target `cor` dimnames must be exactly the fitted ",
               "effect modifiers.", call. = FALSE)
        }
        cor <- cor[ems, ems, drop = FALSE]
      }
      if (!is.numeric(cor) || !identical(dim(cor), c(length(ems), length(ems))) ||
          any(!is.finite(cor)) || !isTRUE(all.equal(cor, t(cor))) ||
          any(abs(diag(cor) - 1) > 1e-8) ||
          inherits(try(chol(cor), silent = TRUE), "try-error")) {
        stop("Summary target `cor` must be a finite positive-definite ",
             "correlation matrix with one row and column per modifier.",
             call. = FALSE)
      }
    }
    resolved <- .cpaic_marginal_summary_nodes(
      object, means, sds, tmargins, cor, n_int_integer
    )
    X <- resolved$X
    w <- resolved$weights
    provenance <- list(
      type = "summary",
      specification = list(means = means, sds = sds, margins = tmargins,
                           cor = cor, n_int = n_int_integer),
      integration = resolved$integration,
      realized_means = resolved$realized_means %||% means,
      realized_sds = resolved$realized_sds %||% sds,
      realized_cor = resolved$realized_cor,
      reference_cor = resolved$reference_cor
    )
  } else {
    if (!is.data.frame(target)) {
      stop("`target` must be a data frame or a summary-distribution list.",
           call. = FALSE)
    }
    if (nrow(target) < 1L) {
      stop("An empirical `target` must contain at least one row.",
           call. = FALSE)
    }
    miss <- setdiff(ems, names(target))
    if (length(miss)) {
      stop("`target` is missing effect modifier(s): ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    }
    if (!setequal(names(target), ems) || length(names(target)) != length(ems)) {
      stop("An empirical `target` must contain exactly the fitted effect ",
           "modifiers: ", paste(ems, collapse = ", "), ".", call. = FALSE)
    }
    bad_type <- ems[!vapply(target[ems], function(x)
      is.numeric(x) && !is.factor(x), logical(1))]
    if (length(bad_type)) {
      stop("Target effect modifiers must be numeric: ",
           paste(bad_type, collapse = ", "), ".", call. = FALSE)
    }
    X <- as.matrix(target[, ems, drop = FALSE])
    storage.mode(X) <- "double"
    if (any(!is.finite(X))) {
      stop("Target effect-modifier values must be finite.", call. = FALSE)
    }
    for (j in seq_along(ems)) {
      v <- ems[j]
      if (margins[[v]] == "bernoulli" && any(!X[, j] %in% c(0, 1))) {
        stop("Empirical Bernoulli target values for `", v,
             "` must be actual 0 or 1 rows.", call. = FALSE)
      }
      if (margins[[v]] %in% c("gamma", "lognormal") && any(X[, j] <= 0)) {
        stop("Empirical target values for `", v, "` must be positive.",
             call. = FALSE)
      }
      if (margins[[v]] == "beta" && any(X[, j] <= 0 | X[, j] >= 1)) {
        stop("Empirical beta target values for `", v, "` must lie in (0, 1).",
             call. = FALSE)
      }
    }
    if (is.null(weights)) weights <- rep(1, nrow(X))
    if (!is.numeric(weights) || length(weights) != nrow(X) ||
        any(!is.finite(weights)) || any(weights < 0) ||
        !any(weights > 0)) {
      stop("`weights` must have one finite nonnegative value per target row ",
           "and a positive total.", call. = FALSE)
    }
    scaled_weights <- as.numeric(weights) / max(weights)
    w <- scaled_weights / sum(scaled_weights)
    provenance <- list(type = "empirical", specification = NULL,
                       integration = "weighted empirical rows",
                       realized_cor = NULL, reference_cor = NULL)
  }

  keep <- w > 0
  X <- X[keep, , drop = FALSE]
  w <- w[keep] / sum(w[keep])
  list(X = X, weights = w, provenance = provenance)
}

.cpaic_marginal_draw_blocks <- function(object) {
  C <- object$C.matrix
  Q <- length(object$effect_modifiers)
  draws <- .cpaic_draws_matrix(
    object$fit, c("beta", "breg", "gamma", "mu")
  )
  bcols <- paste0("beta[", seq_len(ncol(C)), "]")
  rcols <- paste0("breg[", seq_len(Q), "]")
  mcols <- if (length(object$study_levels)) {
    paste0("mu[", seq_along(object$study_levels), "]")
  } else {
    character()
  }
  if (is.null(colnames(draws)) ||
      !all(c(bcols, rcols, mcols) %in% colnames(draws)) ||
      nrow(draws) < 1L) {
    stop("The fit does not contain compatible marginal-prediction draws.",
         call. = FALSE)
  }
  B <- draws[, bcols, drop = FALSE]
  R <- draws[, rcols, drop = FALSE]
  M <- draws[, mcols, drop = FALSE]
  G <- array(NA_real_, dim = c(nrow(B), ncol(C), Q))
  for (q in seq_len(Q)) {
    cols <- paste0("gamma[", seq_len(ncol(C)), ",", q, "]")
    if (!all(cols %in% colnames(draws))) {
      stop("The fit does not contain all component x modifier draws.",
           call. = FALSE)
    }
    G[, , q] <- draws[, cols, drop = FALSE]
  }
  if (any(!is.finite(B)) || any(!is.finite(R)) || any(!is.finite(G)) ||
      any(!is.finite(M))) {
    stop("Marginal prediction requires finite posterior parameter draws.",
         call. = FALSE)
  }
  list(beta = B, breg = R, gamma = G, mu = M)
}

.cpaic_marginal_mu <- function(
    object, baseline_study, required = TRUE, draws = NULL) {
  if (!required && is.null(baseline_study)) return(NULL)
  studies <- object$study_levels
  if (is.null(studies) || !length(studies)) {
    stop("The fit has no retained study-level metadata needed to select a ",
         "baseline study.", call. = FALSE)
  }
  if (!is.character(baseline_study) || length(baseline_study) != 1L ||
      is.na(baseline_study) || !baseline_study %in% studies) {
    stop("`baseline_study` must name one of: ", paste(studies, collapse = ", "),
         ".", call. = FALSE)
  }
  M <- draws %||% .cpaic_draws_matrix(object$fit, "mu")
  idx <- match(baseline_study, studies)
  col <- paste0("mu[", idx, "]")
  if (!col %in% colnames(M)) {
    stop("The fit does not contain the selected study intercept draws.",
         call. = FALSE)
  }
  out <- M[, col]
  if (!length(out) || any(!is.finite(out))) {
    stop("The selected study intercept draws are missing or non-finite.",
         call. = FALSE)
  }
  out
}

.cpaic_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.cpaic_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

.cpaic_marginal_predict_non_survival <- function(
    object, X, weights, treatments, baseline_study, measure) {
  blocks <- .cpaic_marginal_draw_blocks(object)
  C <- object$C.matrix[treatments, , drop = FALSE]
  nd <- nrow(blocks$beta)
  nn <- nrow(X)
  nt <- length(treatments)
  base_theta <- blocks$beta %*% t(C)
  prognostic <- blocks$breg %*% t(X)
  logw <- log(weights)

  theta <- array(0, dim = c(nd, nn, nt))
  for (t in seq_len(nt)) {
    theta[, , t] <- matrix(base_theta[, t], nd, nn)
    for (q in seq_len(ncol(X))) {
      Gq <- matrix(blocks$gamma[, , q], nrow = nd, ncol = ncol(C))
      tq <- Gq %*% C[t, ]
      theta[, , t] <- theta[, , t] + tcrossprod(tq, X[, q])
    }
  }

  if (object$family == "gaussian") {
    values <- matrix(NA_real_, nd, nt, dimnames = list(NULL, treatments))
    for (t in seq_len(nt)) {
      theta_t <- matrix(theta[, , t], nrow = nd, ncol = nn)
      values[, t] <- as.numeric(theta_t %*% weights)
    }
    return(list(values = values, log_values = NULL, log_complements = NULL))
  }

  need_mu <- object$family == "binomial" || identical(measure, "rate_difference")
  mu <- .cpaic_marginal_mu(
    object, baseline_study, required = need_mu, draws = blocks$mu
  )
  if (is.null(mu)) mu <- rep(0, nd)
  if (length(mu) != nd) {
    stop("The study intercept and treatment-effect blocks have inconsistent ",
         "draw counts.", call. = FALSE)
  }

  if (object$family == "poisson") {
    log_values <- matrix(NA_real_, nd, nt, dimnames = list(NULL, treatments))
    for (t in seq_len(nt)) {
      eta <- prognostic + theta[, , t, drop = FALSE][, , 1L] + mu
      for (d in seq_len(nd)) {
        log_values[d, t] <- .cpaic_logsumexp(logw + eta[d, ])
      }
    }
    return(list(values = NULL, log_values = log_values,
                log_complements = NULL))
  }

  log_values <- log_complements <- matrix(
    NA_real_, nd, nt, dimnames = list(NULL, treatments))
  for (t in seq_len(nt)) {
    eta <- prognostic + theta[, , t, drop = FALSE][, , 1L] + mu
    for (d in seq_len(nd)) {
      log_values[d, t] <- .cpaic_logsumexp(
        logw - .cpaic_softplus(-eta[d, ]))
      log_complements[d, t] <- .cpaic_logsumexp(
        logw - .cpaic_softplus(eta[d, ]))
    }
  }
  list(values = NULL, log_values = log_values,
       log_complements = log_complements)
}

.cpaic_marginal_basis <- function(object, t1, t2, X, weights) {
  C <- object$C.matrix
  Njoint <- .cpaic_null_space(object$joint_design)
  Nipd <- .cpaic_null_space(object$joint_design_ipd)
  if (object$family == "gaussian") {
    xbar <- as.numeric(crossprod(weights, X))
    v <- .cpaic_target_vec(C[t1, ] - C[t2, ], xbar)
    if (!.cpaic_in_rowspace(matrix(v, 1L), Njoint)) return("not identified")
    if (.cpaic_in_rowspace(matrix(v, 1L), Nipd)) return("exact")
    return("first-order screen")
  }
  ok_treatment <- function(trt) {
    V <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
      .cpaic_target_vec(C[trt, ], X[i, ])
    }))
    all(.cpaic_in_rowspace(V, Njoint))
  }
  if (!ok_treatment(t1) || !ok_treatment(t2))
    "first-order screen failed" else
    "first-order screen"
}

.cpaic_marginal_pairs <- function(trts, reference, all_contrasts) {
  if (all_contrasts) {
    p <- expand.grid(t1 = trts, t2 = trts, stringsAsFactors = FALSE)
    p[p$t1 != p$t2, , drop = FALSE]
  } else {
    data.frame(t1 = setdiff(trts, reference), t2 = reference,
               stringsAsFactors = FALSE)
  }
}

.cpaic_marginal_signed_exp_difference <- function(log_x, log_y) {
  out <- rep(NA_real_, length(log_x))
  equal <- log_x == log_y
  out[equal] <- 0
  x_zero <- is.infinite(log_x) & log_x < 0 & is.finite(log_y)
  y_zero <- is.infinite(log_y) & log_y < 0 & is.finite(log_x)
  out[x_zero] <- -exp(log_y[x_zero])
  out[y_zero] <- exp(log_x[y_zero])
  unequal <- !equal & is.finite(log_x) & is.finite(log_y)
  if (any(unequal)) {
    x_larger <- log_x[unequal] > log_y[unequal]
    larger <- pmax(log_x[unequal], log_y[unequal])
    smaller <- pmin(log_x[unequal], log_y[unequal])
    log_absolute <- larger + log(-expm1(smaller - larger))
    magnitude <- ifelse(
      log_absolute > log(.Machine$double.xmax), Inf, exp(log_absolute)
    )
    out[unequal] <- ifelse(x_larger, magnitude, -magnitude)
  }
  out
}

.cpaic_marginal_contrast <- function(pred, t1, t2, measure) {
  ratio <- grepl("(_ratio|odds_ratio)$", measure) ||
    identical(measure, "time_specific_hazard_ratio")
  if (!is.null(pred$log_values)) {
    l1 <- pred$log_values[, t1]
    l2 <- pred$log_values[, t2]
    if (measure == "survival_ratio" &&
        !is.null(pred$log_neg_log_values)) {
      out <- l1 - l2
      unrepresentable <- !is.finite(out)
      if (any(unrepresentable)) {
        a1 <- pred$log_neg_log_values[, t1]
        a2 <- pred$log_neg_log_values[, t2]
        out[unrepresentable] <- .cpaic_marginal_signed_exp_difference(
          a2[unrepresentable], a1[unrepresentable]
        )
      }
      return(out)
    }
    if (measure == "odds_ratio") {
      return((l1 - pred$log_complements[, t1]) -
               (l2 - pred$log_complements[, t2]))
    }
    if (measure == "risk_difference" &&
        !is.null(pred$log_complements)) {
      return(.cpaic_marginal_signed_exp_difference(
        pred$log_complements[, t2], pred$log_complements[, t1]
      ))
    }
    if (ratio) return(l1 - l2)
    return(.cpaic_marginal_signed_exp_difference(l1, l2))
  }
  v1 <- pred$values[, t1]
  v2 <- pred$values[, t2]
  if (ratio) log(v1) - log(v2) else v1 - v2
}

.cpaic_marginal_summary_row <- function(
    draws, t1, t2, measure, level, backtransf, basis, time = NULL) {
  ratio <- grepl("(_ratio|odds_ratio)$", measure) ||
    identical(measure, "time_specific_hazard_ratio")
  if (!length(draws) || any(!is.finite(draws))) {
    stop("Marginal contrast draws are missing or non-finite.", call. = FALSE)
  }
  natural <- if (ratio) exp(draws) else draws
  reported <- if (ratio && !backtransf) draws else natural
  a <- (1 - level) / 2
  if (basis %in% c("not identified", "first-order screen failed")) {
    reported[] <- NA_real_
    draws[] <- NA_real_
  }
  q <- if (all(is.na(reported))) c(NA_real_, NA_real_) else
    stats::quantile(reported, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  out <- data.frame(
    treatment = t1, comparator = t2,
    estimate = mean(reported, na.rm = TRUE),
    estimate_contrast = mean(draws, na.rm = TRUE),
    se_contrast = stats::sd(draws, na.rm = TRUE),
    contrast_scale = if (ratio) "log ratio" else "natural difference",
    lower = q[1], upper = q[2],
    scale = if (ratio && !backtransf) "link" else "natural",
    pr_gt0 = mean(draws > 0, na.rm = TRUE), basis = basis,
    stringsAsFactors = FALSE)
  if (all(is.na(draws))) {
    out[c("estimate", "estimate_contrast", "se_contrast", "pr_gt0")] <-
      NA_real_
  }
  if (!is.null(time)) out$time <- time
  out
}

#' @export
marginal_effects.cpaic_mlnmr <- function(
    object, target, weights = NULL, reference = NULL,
    all_contrasts = FALSE, measure = NULL, baseline_study = NULL,
    times = NULL, random_effect = "population", backtransf = TRUE,
    level = 0.95, ...) {
  dots <- list(...)
  if (length(dots)) {
    nms <- names(dots)
    label <- if (is.null(nms) || any(!nzchar(nms))) "unnamed arguments" else
      paste0("`", paste(nms, collapse = "`, `"), "`")
    stop("Unused argument(s): ", label, ".", call. = FALSE)
  }
  all_contrasts <- .cpaic_marginal_scalar_flag(all_contrasts, "all_contrasts")
  backtransf <- .cpaic_marginal_scalar_flag(backtransf, "backtransf")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!identical(random_effect, "population")) {
    stop("`random_effect` currently supports only \"population\", which sets ",
         "study-arm deviations to zero.", call. = FALSE)
  }
  trts <- rownames(object$C.matrix)
  if (is.null(reference)) reference <- object$reference
  .cpaic_check_ref_level(reference, trts, level)

  allowed <- switch(
    object$family,
    binomial = c("odds_ratio", "risk_ratio", "risk_difference"),
    gaussian = "mean_difference",
    poisson = c("rate_ratio", "rate_difference"),
    survival = c("survival_difference", "survival_ratio",
                 "risk_difference", "risk_ratio", "rmst_difference",
                 "rmst_ratio", "time_specific_hazard_ratio"),
    stop("Unsupported cML-NMR family.", call. = FALSE))
  if (is.null(measure)) {
    if (object$family == "survival") {
      stop("Survival marginal effects require an explicit `measure`.",
           call. = FALSE)
    }
    measure <- allowed[1L]
  }
  if (is.character(measure) && length(measure) == 1L &&
      identical(measure, "hazard_ratio")) {
    stop("Use `measure = \"time_specific_hazard_ratio\"` with explicit times.",
         call. = FALSE)
  }
  if (!is.character(measure) || length(measure) != 1L ||
      is.na(measure) || !measure %in% allowed) {
    stop("`measure` for a ", object$family, " fit must be one of: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  target_data <- .cpaic_marginal_target(object, target, weights)
  X <- target_data$X
  w <- target_data$weights
  pairs <- .cpaic_marginal_pairs(trts, reference, all_contrasts)

  if (object$family == "survival") {
    if (is.null(times) || !is.numeric(times) || !length(times) ||
        any(!is.finite(times)) || any(times <= 0)) {
      stop("`times` must contain positive finite survival prediction times.",
           call. = FALSE)
    }
    type <- if (grepl("^survival_", measure)) "survival" else
      if (grepl("^risk_", measure)) "risk" else
        if (grepl("^rmst_", measure)) "rmst" else "hazard"
    pred <- .cpaic_marginal_survival_predictions(
      object, X, w, trts, baseline_study, times, type)
    rows <- list()
    z <- 1L
    for (i in seq_len(nrow(pairs))) {
      basis <- .cpaic_marginal_basis(
        object, pairs$t1[i], pairs$t2[i], X, w)
      for (j in seq_along(pred$index)) {
        values <- matrix(
          pred$values[, , j], nrow = dim(pred$values)[1L],
          ncol = dim(pred$values)[2L]
        )
        log_values <- matrix(
          pred$log_values[, , j], nrow = dim(pred$log_values)[1L],
          ncol = dim(pred$log_values)[2L]
        )
        log_neg_log_values <- if (is.null(pred$log_neg_log_values)) {
          NULL
        } else {
          matrix(
            pred$log_neg_log_values[, , j],
            nrow = dim(pred$log_neg_log_values)[1L],
            ncol = dim(pred$log_neg_log_values)[2L]
          )
        }
        log_complements <- if (is.null(pred$log_complements)) {
          NULL
        } else {
          matrix(
            pred$log_complements[, , j],
            nrow = dim(pred$log_complements)[1L],
            ncol = dim(pred$log_complements)[2L]
          )
        }
        colnames(values) <- trts
        colnames(log_values) <- trts
        if (!is.null(log_neg_log_values)) {
          colnames(log_neg_log_values) <- trts
        }
        if (!is.null(log_complements)) {
          colnames(log_complements) <- trts
        }
        d <- .cpaic_marginal_contrast(
          list(values = values, log_values = log_values,
               log_complements = log_complements,
               log_neg_log_values = log_neg_log_values),
          pairs$t1[i], pairs$t2[i], measure)
        rows[[z]] <- .cpaic_marginal_summary_row(
          d, pairs$t1[i], pairs$t2[i], measure, level, backtransf,
          basis, time = pred$index[j])
        z <- z + 1L
      }
    }
    out <- do.call(rbind, rows)
  } else {
    if (!is.null(times)) {
      stop("`times` is available only for survival marginal measures.",
           call. = FALSE)
    }
    needs_baseline <- object$family == "binomial" ||
      measure == "rate_difference"
    if (!needs_baseline && !is.null(baseline_study)) {
      stop("`baseline_study` is not used for this marginal measure.",
           call. = FALSE)
    }
    pred <- .cpaic_marginal_predict_non_survival(
      object, X, w, trts, baseline_study, measure)
    rows <- lapply(seq_len(nrow(pairs)), function(i) {
      basis <- .cpaic_marginal_basis(
        object, pairs$t1[i], pairs$t2[i], X, w)
      d <- .cpaic_marginal_contrast(
        pred, pairs$t1[i], pairs$t2[i], measure)
      .cpaic_marginal_summary_row(
        d, pairs$t1[i], pairs$t2[i], measure, level, backtransf, basis)
    })
    out <- do.call(rbind, rows)
  }

  rownames(out) <- NULL
  ratio <- grepl("(_ratio|odds_ratio)$", measure) ||
    identical(measure, "time_specific_hazard_ratio")
  attr(out, "sm") <- switch(
    measure, odds_ratio = "OR", risk_ratio = "RR",
    risk_difference = "RD", mean_difference = "MD",
    rate_ratio = "IRR", rate_difference = "IRD",
    survival_ratio = "SR", survival_difference = "SD",
    rmst_ratio = "RMSTR", rmst_difference = "RMSTD",
    time_specific_hazard_ratio = "MHR(t)")
  attr(out, "backtransf") <- ratio && backtransf
  attr(out, "scale") <- if (ratio && !backtransf) "link" else "natural"
  attr(out, "target_mean") <- stats::setNames(
    as.numeric(crossprod(w, X)), object$effect_modifiers)
  attr(out, "target") <- attr(out, "target_mean")
  attr(out, "target_distribution") <- list(
    type = target_data$provenance$type,
    nodes = as.data.frame(X),
    weights = w,
    specification = target_data$provenance$specification,
    integration = target_data$provenance$integration,
    realized_means = target_data$provenance$realized_means %||%
      attr(out, "target_mean"),
    realized_sds = target_data$provenance$realized_sds %||% NULL,
    realized_cor = target_data$provenance$realized_cor %||% NULL,
    reference_cor = target_data$provenance$reference_cor %||% NULL
  )
  attr(out, "contrast_scale") <- if (ratio) "log ratio" else
    "natural difference"
  attr(out, "estimand") <- "marginal"
  attr(out, "measure") <- measure
  attr(out, "baseline_study") <- baseline_study
  attr(out, "random_effect") <- "population"
  attr(out, "estimability_scope") <- "treatment beta/Gamma surface only"
  attr(out, "diagnostic_status") <- object$diagnostics$status %||% "unknown"
  class(out) <- c("cpaic_effects", "data.frame")
  out
}
