# Component-additive ML-NMR (Phase 2, Bayesian) ------------------------------

#' Build numerical integration points for an aggregate study
#'
#' Sobol' quasi-Monte-Carlo points mapped through normal margins defined by
#' the study covariate means and standard deviations. Effect modifiers are
#' assumed continuous in this first version (a Gaussian-copula correlation
#' from IPD is a planned extension).
#' @noRd
.cpaic_integration_points <- function(means, sds, n_int) {
  P <- length(means)
  if (!requireNamespace("randtoolbox", quietly = TRUE)) {
    stop("Package 'randtoolbox' is required for cmlnmr() integration.",
         call. = FALSE)
  }
  u <- randtoolbox::sobol(n = n_int, dim = P)
  u <- matrix(u, nrow = n_int, ncol = P)
  X <- matrix(NA_real_, n_int, P)
  for (p in seq_len(P)) {
    X[, p] <- stats::qnorm(u[, p], mean = means[p], sd = sds[p])
  }
  X
}

#' Component-additive multilevel network meta-regression (ML-NMR)
#'
#' The Bayesian flagship of cpaic. The relative effect of every treatment is
#' the sum of its component effects (`theta = C beta`), estimated jointly
#' from individual patient data and aggregate data. Aggregate arms are
#' fitted by integrating the individual-level model over each study's
#' covariate distribution (numerical integration). Because disconnected
#' sub-networks share component parameters, the network is connected by
#' construction; population-adjusted relative effects in any target
#' population follow from the fitted component effects and effect-modifier
#' interactions.
#'
#' This first version supports binary outcomes with continuous effect
#' modifiers and is fitted with `cmdstanr`. Other families and a
#' Gaussian-copula integration are planned.
#'
#' Note on identifiability: a component whose effect-modifier interaction is
#' informed only by aggregate arms (no individual patient data) is weakly
#' identified, because the main and interaction effects are only constrained
#' through the integrated population-average outcome. Such components rely on
#' the regression prior (`prior_reg_sd`); supply individual patient data for
#' those components where possible.
#'
#' @param ipd Individual patient data: a data frame with a study column, a
#'   treatment column, a binary outcome, and the effect-modifier columns.
#' @param agd Aggregate data, one row per arm, with a study column, a
#'   treatment column, event count `r`, sample size `n`, and, for each
#'   effect modifier `x`, columns `x_mean` and `x_sd`.
#' @param effect_modifiers Character vector of effect-modifier names.
#' @param inactive,sep.comps Component coding (see [cpaic_network()]).
#' @param family Outcome family. Currently `"binomial"`.
#' @param study,trt,outcome,r,n Column names.
#' @param n_int Number of integration points per aggregate arm.
#' @param prior_intercept_sd,prior_beta_sd,prior_reg_sd Prior standard
#'   deviations.
#' @param chains,iter_warmup,iter_sampling,seed Passed to `cmdstanr`.
#' @param ... Reserved.
#'
#' @return An object of class `cpaic_mlnmr` with the `cmdstanr` fit, the
#'   component design, and a tidy table of component effects.
#' @seealso [cmaic()], [cstc()], [cnma_bridge()]
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE)
#' \donttest{
#' # IPD (study, treatment, binary outcome, effect modifier) and aggregate
#' # arms with per-study covariate summaries (x1_mean, x1_sd).
#' ipd <- data.frame(.study = "S1",
#'                   .trt = rep(c("Placebo", "A"), each = 100),
#'                   .y = rbinom(200, 1, 0.5), x1 = rnorm(200))
#' agd <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"),
#'                   r = c(40, 55), n = c(100, 100),
#'                   x1_mean = c(0.2, 0.2), x1_sd = c(1, 1))
#' fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
#'               chains = 2, iter_warmup = 200, iter_sampling = 200)
#' component_effects(fit)
#' }
#' @export
cmlnmr <- function(ipd, agd, effect_modifiers, inactive = NULL,
                   sep.comps = "+", family = "binomial",
                   study = ".study", trt = ".trt", outcome = ".y",
                   r = "r", n = "n", n_int = 64L,
                   prior_intercept_sd = 10, prior_beta_sd = 10,
                   prior_reg_sd = 2.5, chains = 4L, iter_warmup = 500L,
                   iter_sampling = 500L, seed = NULL, ...) {
  family <- match.arg(family, c("binomial"))
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmlnmr() needs 'cmdstanr'. Install it from ",
         "https://stan-dev.r-universe.dev .", call. = FALSE)
  }
  ipd <- as.data.frame(ipd)
  agd <- as.data.frame(agd)
  Q <- length(effect_modifiers)

  # Component design across all treatments.
  all_trts <- sort(unique(c(as.character(ipd[[trt]]),
                            as.character(agd[[trt]]))))
  C <- build_C_matrix(all_trts, sep.comps = sep.comps, inactive = inactive)
  comps <- colnames(C)

  studies <- sort(unique(c(as.character(ipd[[study]]),
                           as.character(agd[[study]]))))
  sidx <- function(s) match(as.character(s), studies)

  # IPD design.
  Tc_ipd <- C[match(as.character(ipd[[trt]]), rownames(C)), , drop = FALSE]
  X_ipd <- as.matrix(ipd[, effect_modifiers, drop = FALSE])
  y_ipd <- as.integer(ipd[[outcome]])

  # AgD integration points (per arm), stacked.
  X_list <- vector("list", nrow(agd))
  for (a in seq_len(nrow(agd))) {
    means <- vapply(effect_modifiers,
                    function(v) agd[[paste0(v, "_mean")]][a], numeric(1))
    sds <- vapply(effect_modifiers,
                  function(v) agd[[paste0(v, "_sd")]][a], numeric(1))
    X_list[[a]] <- .cpaic_integration_points(means, sds, n_int)
  }
  X_agd_int <- do.call(rbind, X_list)
  Tc_agd <- C[match(as.character(agd[[trt]]), rownames(C)), , drop = FALSE]

  standata <- list(
    N_ipd = nrow(ipd), N_agd = nrow(agd), N_studies = length(studies),
    C = ncol(C), P = Q, Q = Q, n_int = as.integer(n_int),
    y_ipd = y_ipd, study_ipd = sidx(ipd[[study]]),
    Tc_ipd = Tc_ipd, X_ipd = X_ipd,
    em_idx = matrix(rep(seq_len(Q), each = nrow(ipd)), nrow = nrow(ipd)),
    r_agd = as.integer(agd[[r]]), n_agd = as.integer(agd[[n]]),
    study_agd = sidx(agd[[study]]), Tc_agd = Tc_agd,
    X_agd_int = X_agd_int,
    prior_intercept_sd = prior_intercept_sd,
    prior_beta_sd = prior_beta_sd, prior_reg_sd = prior_reg_sd
  )

  mod <- .cpaic_stan_model(family)
  fit <- mod$sample(
    data = standata, chains = chains, parallel_chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed %||% 1L, refresh = 0, show_messages = FALSE)

  beta_draws <- fit$draws("beta", format = "draws_matrix")
  comp_tbl <- data.frame(
    component = comps,
    estimate = apply(beta_draws, 2, mean),
    se = apply(beta_draws, 2, stats::sd),
    lower = apply(beta_draws, 2, stats::quantile, 0.025),
    upper = apply(beta_draws, 2, stats::quantile, 0.975),
    row.names = NULL, stringsAsFactors = FALSE
  )

  structure(
    list(fit = fit, components = comp_tbl, C.matrix = C, comps = comps,
         family = family, effect_modifiers = effect_modifiers,
         method = "cML-NMR"),
    class = c("cpaic_mlnmr", "cpaic_fit")
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Compile (and cache) a cpaic Stan model with cmdstanr
#'
#' The model is compiled into a per-user cache directory rather than next to
#' the installed `.stan` file, so the package source/installation tree never
#' acquires an executable. The compiled binary is reused across calls.
#' @noRd
.cpaic_stan_model <- function(family) {
  stan_file <- system.file("stan", paste0("cpaic_", family, ".stan"),
                           package = "cpaic")
  if (stan_file == "") {
    # devtools::load_all() path
    stan_file <- file.path("inst", "stan", paste0("cpaic_", family, ".stan"))
  }
  cache <- tools::R_user_dir("cpaic", "cache")
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)
  dest <- file.path(cache, basename(stan_file))
  if (!file.exists(dest) ||
      file.mtime(stan_file) > file.mtime(dest)) {
    file.copy(stan_file, dest, overwrite = TRUE)
  }
  cmdstanr::cmdstan_model(dest)
}

#' @export
component_effects.cpaic_mlnmr <- function(object, ...) object$components

#' @export
print.cpaic_mlnmr <- function(x, ...) {
  cat("cpaic: component-additive ML-NMR (Bayesian, ", x$family, ")\n",
      sep = "")
  cat("  Effect modifiers: ", paste(x$effect_modifiers, collapse = ", "),
      "\n\n", sep = "")
  comp <- x$components
  comp[, c("estimate", "se", "lower", "upper")] <-
    round(comp[, c("estimate", "se", "lower", "upper")], 3)
  print(comp, row.names = FALSE)
  invisible(x)
}
