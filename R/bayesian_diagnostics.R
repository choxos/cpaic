# Bayesian diagnostics and prior checks

#' Extract pointwise log-likelihood draws
#' @noRd
.cpaic_log_lik <- function(object) {
  if (!inherits(object, "cpaic_mlnmr")) {
    stop("`object` must be a cpaic_mlnmr fit.", call. = FALSE)
  }
  out <- .cpaic_draws_matrix(object$fit, "log_lik")
  if (!ncol(out)) {
    stop("The fit does not contain pointwise `log_lik` draws.", call. = FALSE)
  }
  out
}

#' Pareto-smoothed importance sampling leave-one-out cross-validation
#'
#' This is **observation-level** LOO: it measures within-study interpolation
#' (leaving out one IPD patient or one reconstructed pseudo-observation), which
#' is not the scientific question in a disconnected network. It does **not**
#' validate cross-gap prediction (a new study, a held-out treatment contrast, or
#' a held-out sub-network); a good pointwise LOO can coexist with a wrong
#' cross-gap extrapolation. Grouped leave-one-study-out is not yet implemented,
#' so `unit` accepts only `"observation"`.
#'
#' @param x A [cmlnmr()] fit.
#' @param unit Predictive unit; only `"observation"` is supported.
#' @param ... Passed to [loo::loo.matrix()].
#' @return A `psis_loo` object from the `loo` package.
#' @importFrom loo loo
#' @export
loo.cpaic_mlnmr <- function(x, unit = "observation", ...) {
  unit <- match.arg(unit, "observation")
  loo::loo(.cpaic_log_lik(x), ...)
}

#' Widely applicable information criterion
#'
#' @param x A [cmlnmr()] fit.
#' @param ... Passed to [loo::waic.matrix()].
#' @return A `waic` object from the `loo` package.
#' @importFrom loo waic
#' @export
waic.cpaic_mlnmr <- function(x, ...) {
  loo::waic(.cpaic_log_lik(x), ...)
}

#' Deviance information criterion
#'
#' Computes DIC with the variance penalty `pV`, following the survival-model
#' implementation in multinma. The pointwise deviance is `-2 * log_lik`; the
#' effective parameter count is half the posterior variance of total deviance.
#'
#' @param x A fitted model.
#' @param ... Unused.
#' @return For [cmlnmr()] fits, a `cpaic_dic` object with DIC, mean deviance,
#'   the `pV` penalty, and pointwise mean deviance.
#' @export
dic <- function(x, ...) {
  UseMethod("dic")
}

#' @export
dic.cpaic_mlnmr <- function(x, ...) {
  log_lik <- .cpaic_log_lik(x)
  deviance_draws <- -2 * rowSums(log_lik)
  mean_deviance <- mean(deviance_draws)
  pV <- stats::var(deviance_draws) / 2
  structure(
    list(
      dic = mean_deviance + pV,
      penalty = "pV",
      p_eff = pV,
      mean_deviance = mean_deviance,
      pointwise = -2 * colMeans(log_lik)
    ),
    class = "cpaic_dic"
  )
}

#' @export
print.cpaic_dic <- function(x, digits = 1, ...) {
  cat("Deviance information criterion\n")
  cat("  DIC: ", round(x$dic, digits), "\n", sep = "")
  cat("  Mean deviance: ", round(x$mean_deviance, digits), "\n", sep = "")
  cat("  Effective parameters (pV): ", round(x$p_eff, digits), "\n",
      sep = "")
  invisible(x)
}

#' Summarize Poisson replication overflow counts
#'
#' The Poisson Stan model records how many generated counts could not be drawn
#' safely in each iteration. This helper validates those generated quantities
#' before any replicated outcome is summarized.
#' @noRd
.cpaic_poisson_overflow_summary <- function(fit, variables) {
  sources <- c("ipd", "agd")
  valid_variables <- is.list(variables) &&
    all(sources %in% names(variables)) &&
    all(vapply(variables[sources], function(x) {
      is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }, logical(1)))
  if (!valid_variables) {
    stop("Poisson replication overflow metadata are missing or malformed. ",
         "Refit with the current cpaic version before using replicated ",
         "outcomes.", call. = FALSE)
  }

  rows <- lapply(sources, function(source) {
    variable <- variables[[source]]
    counts <- tryCatch(
      .cpaic_draws_matrix(fit, variable),
      error = function(e) {
        stop("Could not read Poisson replication overflow counts from `",
             variable, "`: ", conditionMessage(e),
             ". Refit with the current cpaic version before using replicated ",
             "outcomes.", call. = FALSE)
      }
    )
    if (ncol(counts) != 1L) {
      stop("Poisson replication overflow variable `", variable,
           "` must contain one count per draw.", call. = FALSE)
    }
    counts <- as.numeric(counts)
    if (!length(counts) || any(!is.finite(counts)) || any(counts < 0) ||
        any(counts != floor(counts))) {
      stop("Poisson replication overflow variable `", variable,
           "` contains invalid counts.", call. = FALSE)
    }
    data.frame(
      source = source,
      draws = length(counts),
      draws_with_overflow = sum(counts > 0),
      overflowed_values = sum(counts),
      max_overflow_per_draw = max(counts),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Summarize a prior-predictive cML-NMR fit
#'
#' `cmlnmr(prior_predictive = TRUE)` samples from the prior without adding the
#' observed likelihood. This helper compares a simple statistic of the
#' observed outcomes with the corresponding replicated outcomes. Survival
#' replications are event-by-observed-time indicators because the censoring
#' process is not modeled. For Poisson models, the check stops if any generated
#' count exceeded Stan's safe random-number range. Affected values are explicit
#' sentinels and are never summarized as data.
#'
#' @param object A [cmlnmr()] fit created with `prior_predictive = TRUE`.
#' @param statistic Either `"mean"` or `"sd"`.
#' @param level Central prior-predictive interval level.
#' @return A data frame with observed and replicated summaries for IPD and AgD.
#' @export
prior_predictive_check <- function(object, statistic = c("mean", "sd"),
                                   level = 0.95) {
  statistic <- match.arg(statistic)
  if (!inherits(object, "cpaic_mlnmr") ||
      !isTRUE(object$prior_predictive)) {
    stop("Use a cmlnmr() fit created with `prior_predictive = TRUE`.",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  stat <- switch(statistic, mean = base::mean, sd = stats::sd)
  alpha <- (1 - level) / 2

  if (identical(object$family, "poisson")) {
    overflow <- .cpaic_poisson_overflow_summary(
      object$fit, object$rep_overflow_variables
    )
    affected <- overflow$overflowed_values > 0
    if (any(affected)) {
      detail <- paste0(
        overflow$source[affected], "=",
        overflow$overflowed_values[affected], collapse = ", "
      )
      stop("Poisson replicated outcomes contain RNG overflow sentinels (",
           detail, "). Refit with more appropriate priors or model scaling ",
           "before running prior_predictive_check().", call. = FALSE)
    }
  }

  summarize_source <- function(source) {
    variable <- object$rep_variables[[source]]
    draws <- .cpaic_draws_matrix(object$fit, variable)
    if (identical(object$family, "poisson") && any(draws < 0)) {
      stop("Poisson replicated outcomes for `", source,
           "` contain negative overflow sentinels even though the recorded ",
           "overflow count is zero. The generated quantities are inconsistent ",
           "and cannot be summarized; refit the model.", call. = FALSE)
    }
    replicated <- apply(draws, 1, stat)
    observed <- stat(object$observed[[source]])
    interval <- stats::quantile(replicated, c(alpha, 0.5, 1 - alpha),
                                names = FALSE)
    data.frame(
      source = source,
      statistic = statistic,
      observed = observed,
      rep_lower = interval[1],
      rep_median = interval[2],
      rep_upper = interval[3],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, lapply(c("ipd", "agd"), summarize_source))
  class(out) <- c("cpaic_prior_predictive", "data.frame")
  out
}

#' Refit cML-NMR under tighter and looser priors
#'
#' Prior movement is an empirical identification diagnostic. Contrasts that
#' move substantially when a weakly identified prior is changed should not be
#' interpreted as data-driven.
#'
#' @param object A [cmlnmr()] fit.
#' @param newdata One row of target effect-modifier means, as for
#'   [relative_effects()].
#' @param prior Which scales to vary: the interaction prior, component-effect
#'   prior, or all configurable scale priors.
#' @param tighter,looser Positive multipliers for the fitted prior scales.
#' @param reference Reference treatment. Defaults to the fit reference.
#' @param ... Named arguments overriding the stored refit call, such as fewer
#'   sampling iterations for a screening run.
#' @return A `cpaic_prior_sensitivity` object containing the movement table and
#'   the tighter and looser fits.
#' @export
prior_sensitivity <- function(object, newdata,
                              prior = c("gamma", "beta", "all"),
                              tighter = 0.5, looser = 2,
                              reference = NULL, ...) {
  if (!inherits(object, "cpaic_mlnmr") || is.null(object$refit_args)) {
    stop("`object` must be a refittable cmlnmr() fit.", call. = FALSE)
  }
  if (isTRUE(attr(object, "redacted"))) {
    stop("`object` has been redacted (raw data removed by redact_fit()); ",
         "prior sensitivity needs the data to refit.", call. = FALSE)
  }
  prior <- match.arg(prior)
  if (any(!is.finite(c(tighter, looser))) || tighter <= 0 || looser <= 0 ||
      tighter >= 1 || looser <= 1) {
    stop("`tighter` must be in (0, 1) and `looser` must exceed 1.",
         call. = FALSE)
  }
  if (is.null(reference)) reference <- object$reference
  overrides <- list(...)
  if (length(overrides) &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))))) {
    stop("Prior-sensitivity refit overrides in `...` must be named.",
         call. = FALSE)
  }

  scale_names <- switch(
    prior,
    gamma = "prior_gamma_scale",
    beta = "prior_beta_sd",
    all = c("prior_intercept_sd", "prior_beta_sd", "prior_sigma_sd",
            "prior_reg_sd", "prior_aux_sd", "prior_gamma_scale",
            "prior_tau_scale")
  )
  # A scale that is not part of the fitted call (e.g. absent from an older
  # refit) is skipped rather than multiplied into an empty value.
  scale_names <- scale_names[
    vapply(scale_names, function(nm) is.numeric(object$refit_args[[nm]]) &&
             length(object$refit_args[[nm]]) == 1L, logical(1))]
  if (!length(scale_names)) {
    stop("No prior scales are available to vary for prior = \"", prior,
         "\"; the movement would be spuriously zero.", call. = FALSE)
  }
  scaled_args <- function(multiplier) {
    args <- object$refit_args
    for (name in scale_names) args[[name]] <- args[[name]] * multiplier
    args$prior_predictive <- FALSE
    if (length(overrides)) args[names(overrides)] <- overrides
    args
  }
  tight_fit <- do.call(cmlnmr, scaled_args(tighter))
  loose_fit <- do.call(cmlnmr, scaled_args(looser))

  raw_contrasts <- function(fit) {
    x <- .cpaic_target_x(newdata, fit$effect_modifiers, fit$margins)
    effects <- .cpaic_beta_at(fit, x) %*% t(fit$C.matrix)
    colnames(effects) <- rownames(fit$C.matrix)
    treatments <- setdiff(colnames(effects), reference)
    stats::setNames(
      vapply(treatments, function(trt) {
        mean(effects[, trt] - effects[, reference])
      }, numeric(1)),
      treatments
    )
  }
  base <- raw_contrasts(object)
  tight <- raw_contrasts(tight_fit)
  loose <- raw_contrasts(loose_fit)
  estimable <- estimable_effects_at(object, newdata = newdata,
                                    reference = reference)
  movement <- data.frame(
    treatment = names(base),
    comparator = reference,
    estimate = unname(base),
    tighter = unname(tight[names(base)]),
    looser = unname(loose[names(base)]),
    move_tighter = unname(abs(tight[names(base)] - base)),
    move_looser = unname(abs(loose[names(base)] - base)),
    stringsAsFactors = FALSE
  )
  movement$max_movement <- pmax(movement$move_tighter,
                                movement$move_looser)
  movement$estimable <- estimable$estimable[
    match(movement$treatment, estimable$treatment)]
  target_mean <- .cpaic_target_x(newdata, object$effect_modifiers,
                                 object$margins)

  structure(
    list(movement = movement, fits = list(tighter = tight_fit,
                                         looser = loose_fit),
         prior = prior, multipliers = c(tighter = tighter, looser = looser),
         target = target_mean, target_mean = target_mean,
         estimand = "average_conditional_link"),
    class = "cpaic_prior_sensitivity"
  )
}

#' @export
print.cpaic_prior_sensitivity <- function(x, digits = 3, ...) {
  cat("cML-NMR prior sensitivity: ", x$prior, " prior\n", sep = "")
  out <- x$movement
  numeric <- vapply(out, is.numeric, logical(1))
  out[numeric] <- lapply(out[numeric], round, digits = digits)
  print(out, row.names = FALSE)
  invisible(x)
}
