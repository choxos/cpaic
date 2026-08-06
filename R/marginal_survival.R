# Marginal survival prediction for cML-NMR -----------------------------------

#' Read and order one posterior block used by survival standardization
#' @noRd
.cpaic_marginal_surv_block <- function(draws, variable, columns) {
  if (is.null(colnames(draws)) || !all(columns %in% colnames(draws))) {
    stop("Posterior block `", variable, "` does not contain the expected ",
         "Stan parameter names.", call. = FALSE)
  }
  out <- draws[, columns, drop = FALSE]
  if (nrow(out) < 1L || any(!is.finite(out))) {
    stop("Posterior block `", variable, "` is empty or non-finite.",
         call. = FALSE)
  }
  out
}

#' Row-wise log-sum-exp
#' @noRd
.cpaic_marginal_surv_log_sum_exp <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  xmax <- apply(x, 1L, max)
  out <- rep(-Inf, nrow(x))
  finite <- is.finite(xmax)
  if (any(finite)) {
    centered <- x[finite, , drop = FALSE] - xmax[finite]
    out[finite] <- xmax[finite] + log(rowSums(exp(centered)))
  }
  out
}

#' Element-wise log(exp(a) + exp(b))
#' @noRd
.cpaic_marginal_surv_log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  both_negative_infinite <- is.infinite(a) & a < 0 &
    is.infinite(b) & b < 0
  out[both_negative_infinite] <- -Inf
  out
}

#' log(1 - exp(-exp(log_cumulative))) without avoidable underflow
#' @noRd
.cpaic_marginal_log_risk <- function(log_cumulative) {
  out <- log_cumulative
  regular <- log_cumulative > -20
  if (any(regular)) {
    cumulative <- exp(log_cumulative[regular])
    out[regular] <- log(-expm1(-cumulative))
  }
  out[is.infinite(log_cumulative) & log_cumulative > 0] <- 0
  out
}

#' Cumulative-hazard differences after cancelling the row minimum
#' @noRd
.cpaic_marginal_relative_cumulative <- function(log_cumulative) {
  row_minimum <- apply(log_cumulative, 1L, min)
  delta <- sweep(log_cumulative, 1L, row_minimum, "-")
  out <- matrix(Inf, nrow = nrow(delta), ncol = ncol(delta))
  at_minimum <- delta == 0
  out[at_minimum] <- 0

  positive <- !at_minimum & is.finite(delta)
  if (any(positive)) {
    log_expm1 <- delta[positive]
    moderate <- log_expm1 < 50
    log_expm1[moderate] <- log(expm1(log_expm1[moderate]))
    log_difference <-
      matrix(row_minimum, nrow = nrow(delta), ncol = ncol(delta))[positive] +
      log_expm1
    representable <- log_difference <= log(.Machine$double.xmax)
    values <- rep(Inf, length(log_difference))
    values[representable] <- exp(log_difference[representable])
    out[positive] <- values
  }
  out
}

#' Construct fixed, population-average linear predictors
#' @noRd
.cpaic_marginal_surv_eta <- function(
    object, X, treatments, study_index, draws) {
  C <- object$C.matrix
  n_comp <- ncol(C)
  n_em <- ncol(X)
  beta <- .cpaic_marginal_surv_block(
    draws, "beta", paste0("beta[", seq_len(n_comp), "]")
  )
  breg <- .cpaic_marginal_surv_block(
    draws, "breg", paste0("breg[", seq_len(n_em), "]")
  )
  mu <- .cpaic_marginal_surv_block(
    draws, "mu", paste0("mu[", seq_along(object$study_levels), "]")
  )
  expected_gamma <- unlist(lapply(seq_len(n_em), function(q) {
    paste0("gamma[", seq_len(n_comp), ",", q, "]")
  }), use.names = FALSE)
  gamma <- .cpaic_marginal_surv_block(
    draws, "gamma", expected_gamma
  )

  n_draw <- nrow(beta)
  if (nrow(breg) != n_draw || nrow(mu) != n_draw || nrow(gamma) != n_draw) {
    stop("Posterior parameter blocks have inconsistent draw counts.",
         call. = FALSE)
  }
  common <- outer(mu[, study_index], rep(1, nrow(X))) +
    breg %*% t(X)

  eta <- lapply(treatments, function(treatment) {
    tc <- C[treatment, , drop = TRUE]
    interaction <- matrix(0, nrow = n_draw, ncol = n_em)
    for (q in seq_len(n_em)) {
      cols <- paste0("gamma[", seq_len(n_comp), ",", q, "]")
      interaction[, q] <- as.numeric(
        gamma[, cols, drop = FALSE] %*% tc
      )
    }
    common + outer(as.numeric(beta %*% tc), rep(1, nrow(X))) +
      interaction %*% t(X)
  })
  names(eta) <- treatments
  eta
}

#' Standardized survival, risk, or hazard at specified times
#' @noRd
.cpaic_marginal_surv_at_times <- function(spec, times, coefficients, eta,
                                           log_weights, type) {
  basis <- .cpaic_survival_basis_eval(
    spec, time = times, start_time = rep(0, length(times)),
    entry_time = rep(0, length(times))
  )
  cumulative_base <- coefficients %*% t(basis$itime)
  hazard_base <- coefficients %*% t(basis$time)
  n_draw <- nrow(coefficients)
  treatments <- names(eta)
  values <- array(
    NA_real_, dim = c(n_draw, length(treatments), length(times)),
    dimnames = list(NULL, treatments, as.character(times))
  )
  log_values <- values
  log_neg_log_values <- if (type == "survival") values else NULL
  log_complements <- if (type == "risk") values else NULL

  for (k in seq_along(times)) {
    log_H0 <- log(cumulative_base[, k])
    log_h0 <- log(hazard_base[, k])
    for (j in seq_along(treatments)) {
      eta_j <- eta[[j]]
      log_cumulative <- outer(log_H0, rep(1, ncol(eta_j))) + eta_j
      cumulative <- exp(log_cumulative)
      log_survival <- .cpaic_marginal_surv_log_sum_exp(
        sweep(-cumulative, 2L, log_weights, "+")
      )

      if (type == "survival") {
        log_value <- log_survival
        relative_cumulative <-
          .cpaic_marginal_relative_cumulative(log_cumulative)
        log_relative_survival <- .cpaic_marginal_surv_log_sum_exp(
          sweep(-relative_cumulative, 2L, log_weights, "+")
        )
        log_minimum_cumulative <- apply(log_cumulative, 1L, min)
        log_correction <- log(-log_relative_survival)
        log_neg_log_values[, j, k] <- .cpaic_marginal_surv_log_add_exp(
          log_minimum_cumulative, log_correction
        )
      } else if (type == "risk") {
        log_individual_risk <- .cpaic_marginal_log_risk(log_cumulative)
        log_value <- .cpaic_marginal_surv_log_sum_exp(
          sweep(log_individual_risk, 2L, log_weights, "+")
        )
        log_complements[, j, k] <- log_survival
      } else {
        log_hazard <- outer(log_h0, rep(1, ncol(eta_j))) + eta_j
        relative_cumulative <-
          .cpaic_marginal_relative_cumulative(log_cumulative)
        log_risk_set <- sweep(-relative_cumulative, 2L, log_weights, "+")
        log_numerator <- .cpaic_marginal_surv_log_sum_exp(
          log_risk_set + log_hazard
        )
        log_denominator <- .cpaic_marginal_surv_log_sum_exp(log_risk_set)
        log_value <- log_numerator - log_denominator
      }
      log_values[, j, k] <- log_value
      values[, j, k] <- exp(log_value)
    }
  }
  list(values = values, log_values = log_values,
       log_neg_log_values = log_neg_log_values,
       log_complements = log_complements)
}

#' Exact RMST under a piecewise-exponential baseline
#' @noRd
.cpaic_marginal_piecewise_rmst <- function(spec, horizons, coefficients, eta,
                                            log_weights) {
  n_draw <- nrow(coefficients)
  treatments <- names(eta)
  values <- array(
    NA_real_, dim = c(n_draw, length(treatments), length(horizons)),
    dimnames = list(NULL, treatments, as.character(horizons))
  )
  log_values <- values

  for (h in seq_along(horizons)) {
    tau <- horizons[h]
    endpoints <- sort(unique(c(0, spec$cut_points[spec$cut_points < tau], tau)))
    lower <- endpoints[-length(endpoints)]
    upper <- endpoints[-1L]
    widths <- upper - lower
    intervals <- findInterval(lower + widths / 2, spec$cut_points) + 1L

    for (j in seq_along(treatments)) {
      eta_j <- eta[[j]]
      cumulative <- matrix(0, nrow = n_draw, ncol = ncol(eta_j))
      log_individual_rmst <- matrix(
        -Inf, nrow = n_draw, ncol = ncol(eta_j)
      )
      for (m in seq_along(widths)) {
        log_rate <- outer(
          log(coefficients[, intervals[m]]), rep(1, ncol(eta_j))
        ) + eta_j
        log_z <- log_rate + log(widths[m])
        z <- exp(log_z)
        log_ratio <- matrix(NA_real_, nrow = nrow(z), ncol = ncol(z))
        small <- log_z < log(1e-5)
        if (any(small)) {
          z_small <- exp(log_z[small])
          log_ratio[small] <- log1p(
            -z_small / 2 + z_small^2 / 6
          )
        }
        large <- log_z > 35
        log_ratio[large] <- -log_z[large]
        moderate <- !small & !large
        if (any(moderate)) {
          z_moderate <- exp(log_z[moderate])
          log_ratio[moderate] <-
            log(-expm1(-z_moderate)) - log_z[moderate]
        }
        log_interval <- -cumulative + log(widths[m]) + log_ratio
        log_individual_rmst <- .cpaic_marginal_surv_log_add_exp(
          log_individual_rmst, log_interval
        )
        cumulative <- cumulative + z
      }
      log_value <- .cpaic_marginal_surv_log_sum_exp(
        sweep(log_individual_rmst, 2L, log_weights, "+")
      )
      log_values[, j, h] <- log_value
      values[, j, h] <- exp(log_value)
    }
  }
  list(values = values, log_values = log_values)
}

#' Treatment-level marginal survival predictions
#'
#' Uses the selected donor study's intercept and baseline hazard. Treatment
#' random effects are set to zero, so predictions describe the population
#' treatment effect rather than a fitted study-arm deviation.
#' @noRd
.cpaic_marginal_survival_predictions <- function(
    object, X, weights, treatments, baseline_study, times,
    type = c("survival", "risk", "rmst", "hazard")) {
  type <- match.arg(type)
  if (!inherits(object, "cpaic_mlnmr") ||
      !identical(object$family, "survival")) {
    stop("Survival marginal predictions require a survival cmlnmr fit.",
         call. = FALSE)
  }
  if (is.null(object$study_levels) || !length(object$study_levels)) {
    stop("The fit does not retain study-level metadata needed for marginal ",
         "survival prediction.", call. = FALSE)
  }
  if (!is.character(baseline_study) || length(baseline_study) != 1L ||
      is.na(baseline_study) ||
      !(baseline_study %in% object$study_levels)) {
    stop("`baseline_study` must name one study retained in the fit.",
         call. = FALSE)
  }
  spec <- object$survival_spec
  if (!inherits(spec, "cpaic_survival_basis") ||
      is.null(spec$max_time) || !is.finite(spec$max_time)) {
    stop("The fit does not retain a valid survival basis specification.",
         call. = FALSE)
  }
  if (type == "rmst" && identical(spec$baseline, "mspline")) {
    stop(
      "Marginal RMST currently requires piecewise-exponential survival ",
      "baselines. M-spline fits still support survival, risk, and ",
      "time-specific marginal hazard measures.",
      call. = FALSE
    )
  }
  study_support <- object$survival_study_support
  donor_support <- if (!is.null(study_support) &&
                       !is.null(names(study_support))) {
    unname(study_support[baseline_study])
  } else {
    NA_real_
  }
  if (length(donor_support) != 1L || !is.finite(donor_support) ||
      donor_support <= 0) {
    stop("The fit does not retain valid follow-up support for the selected ",
         "baseline study.", call. = FALSE)
  }
  times <- as.numeric(times)
  if (!length(times) || any(!is.finite(times)) || any(times <= 0) ||
      any(times > donor_support)) {
    stop("`times` must be positive, finite, and within the fitted survival ",
         "support of baseline study `", baseline_study, "`, ending at ",
         donor_support, ".", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) < 1L ||
      ncol(X) != length(object$effect_modifiers) || any(!is.finite(X))) {
    stop("`X` must be a finite numeric target matrix with one column per ",
         "effect modifier.", call. = FALSE)
  }
  weights <- as.numeric(weights)
  if (length(weights) != nrow(X) || any(!is.finite(weights)) ||
      any(weights <= 0) || !is.finite(sum(weights))) {
    stop("`weights` must contain one positive finite value per target row.",
         call. = FALSE)
  }
  weights <- weights / sum(weights)
  if (!is.character(treatments) || !length(treatments) ||
      anyDuplicated(treatments) ||
      !all(treatments %in% rownames(object$C.matrix))) {
    stop("`treatments` must be unique treatment names in the fitted network.",
         call. = FALSE)
  }

  study_index <- match(baseline_study, object$study_levels)
  draws <- .cpaic_draws_matrix(
    object$fit, c("beta", "breg", "mu", "gamma", "coefficients")
  )
  eta <- .cpaic_marginal_surv_eta(
    object, X, treatments, study_index, draws
  )
  n_base <- spec$n_basis
  coefficient_columns <- paste0(
    "coefficients[", study_index, ",", seq_len(n_base), "]"
  )
  coefficients <- .cpaic_marginal_surv_block(
    draws, "coefficients", coefficient_columns
  )
  if (nrow(coefficients) != nrow(eta[[1L]]) ||
      any(!is.finite(coefficients)) || any(coefficients <= 0)) {
    stop("Survival baseline posterior draws are invalid or have an ",
         "inconsistent draw count.", call. = FALSE)
  }
  log_weights <- log(weights)

  prediction <- if (type == "rmst") {
    .cpaic_marginal_piecewise_rmst(
      spec, times, coefficients, eta, log_weights
    )
  } else {
    .cpaic_marginal_surv_at_times(
      spec, times, coefficients, eta, log_weights, type
    )
  }

  list(values = prediction$values, log_values = prediction$log_values,
       log_neg_log_values = prediction$log_neg_log_values %||% NULL,
       log_complements = prediction$log_complements %||% NULL,
       index = times, outcome = type)
}
