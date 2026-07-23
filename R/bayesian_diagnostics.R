# Bayesian diagnostics and prior checks

#' Extract pointwise log-likelihood draws
#' @noRd
.cpaic_log_lik <- function(object) {
  if (!inherits(object, "cpaic_mlnmr")) {
    stop("`object` must be a cpaic_mlnmr fit.", call. = FALSE)
  }
  out <- as.matrix(object$fit$draws("log_lik", format = "draws_matrix"))
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

#' Summarize a prior-predictive cML-NMR fit
#'
#' `cmlnmr(prior_predictive = TRUE)` samples from the prior without adding the
#' observed likelihood. This helper compares a simple statistic of the
#' observed outcomes with the corresponding replicated outcomes. Survival
#' replications are event-by-observed-time indicators because the censoring
#' process is not modeled.
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

  summarize_source <- function(source) {
    variable <- object$rep_variables[[source]]
    draws <- as.matrix(object$fit$draws(variable, format = "draws_matrix"))
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
#' interpreted as data-driven. This helper reuses the principle in
#' `documentation/validation/estimability_gamma.R`.
#'
#' @param object A [cmlnmr()] fit.
#' @param newdata One target population, as for [relative_effects()].
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
  if (isTRUE(attr(object, "redacted")) || is.null(object$refit_args$ipd)) {
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
    x <- .cpaic_target_x(newdata, fit$effect_modifiers)
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

  structure(
    list(movement = movement, fits = list(tighter = tight_fit,
                                         looser = loose_fit),
         prior = prior, multipliers = c(tighter = tighter, looser = looser),
         target = .cpaic_target_x(newdata, object$effect_modifiers)),
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

