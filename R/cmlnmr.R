# Component-additive ML-NMR (Phase 2, Bayesian) ------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

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
  u <- matrix(randtoolbox::sobol(n = n_int, dim = P), nrow = n_int, ncol = P)
  X <- matrix(NA_real_, n_int, P)
  for (p in seq_len(P)) X[, p] <- stats::qnorm(u[, p], means[p], sds[p])
  X
}

#' Component-additive multilevel network meta-regression (ML-NMR)
#'
#' The Bayesian flagship of cpaic. The relative effect of every treatment is
#' the sum of its component effects (`theta = C beta`), estimated jointly
#' from individual patient data (IPD) and aggregate data (AgD). Aggregate
#' arms are fitted by integrating the individual-level model over each
#' study's covariate distribution. Because disconnected sub-networks share
#' component parameters, the network is connected by construction.
#'
#' Supported families: `"binomial"` (logit), `"gaussian"` (identity),
#' `"poisson"` (log), and `"survival"` (exponential, fitted through the
#' Poisson formulation with person-time exposure). Flexible parametric and
#' spline survival baselines are future work; the frequentist [cstc()] and
#' [cmaic()] cover proportional-hazards survival.
#'
#' Identifiability note: a component whose effect-modifier interaction is
#' informed only by aggregate arms is weakly identified, because main and
#' interaction effects are constrained only through the integrated
#' population-average outcome. Supply IPD for such components where possible;
#' `prior_reg_sd` regularizes otherwise.
#'
#' @param ipd Individual patient data (one row per patient).
#' @param agd Aggregate data (one row per arm) with the per-study covariate
#'   summaries `x_mean`, `x_sd` for each effect modifier `x`.
#' @param effect_modifiers Character vector of effect-modifier names.
#' @param inactive,sep.comps Component coding (see [cpaic_network()]).
#' @param family One of `"binomial"`, `"gaussian"`, `"poisson"`,
#'   `"survival"`.
#' @param study,trt Column names (in both `ipd` and `agd`).
#' @param outcome IPD outcome column: 0/1 (binomial), numeric (gaussian),
#'   count (poisson), or the event indicator for survival.
#' @param time,exposure IPD time / exposure column (survival, poisson).
#' @param r,n,E,se Aggregate columns: events `r`, sample size `n`
#'   (binomial), exposure `E` (poisson/survival), mean `outcome` and its
#'   standard error `se` (gaussian).
#' @param n_int Integration points per aggregate arm (ignored for
#'   `gaussian`, which is exact at the covariate means).
#' @param prior_intercept_sd,prior_beta_sd,prior_reg_sd Prior SDs.
#' @param chains,iter_warmup,iter_sampling,seed Passed to `cmdstanr`.
#' @param ... Reserved.
#'
#' @return An object of class `cpaic_mlnmr` with the `cmdstanr` fit, the
#'   component design, and a tidy table of component effects.
#' @seealso [cmaic()], [cstc()], [cnma_bridge()]
#' @examplesIf requireNamespace("cmdstanr", quietly = TRUE) && !inherits(try(cmdstanr::cmdstan_path(), silent = TRUE), "try-error")
#' \donttest{
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
                   time = ".time", exposure = ".exposure",
                   r = "r", n = "n", E = "E", se = "se", n_int = 64L,
                   prior_intercept_sd = 10, prior_beta_sd = 10,
                   prior_reg_sd = 2.5, chains = 4L, iter_warmup = 500L,
                   iter_sampling = 500L, seed = NULL, ...) {
  family <- match.arg(family,
                      c("binomial", "gaussian", "poisson", "survival"))
  stan_family <- switch(family, binomial = "binomial", gaussian = "normal",
                        poisson = "poisson", survival = "poisson")
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmlnmr() needs 'cmdstanr'. Install from ",
         "https://stan-dev.r-universe.dev .", call. = FALSE)
  }
  ipd <- as.data.frame(ipd)
  agd <- as.data.frame(agd)
  Q <- length(effect_modifiers)

  all_trts <- sort(unique(c(as.character(ipd[[trt]]),
                            as.character(agd[[trt]]))))
  C <- build_C_matrix(all_trts, sep.comps = sep.comps, inactive = inactive)
  comps <- colnames(C)
  studies <- sort(unique(c(as.character(ipd[[study]]),
                           as.character(agd[[study]]))))
  sidx <- function(s) match(as.character(s), studies)

  Tc_ipd <- C[match(as.character(ipd[[trt]]), rownames(C)), , drop = FALSE]
  X_ipd <- as.matrix(ipd[, effect_modifiers, drop = FALSE])
  Tc_agd <- C[match(as.character(agd[[trt]]), rownames(C)), , drop = FALSE]

  # Integration points (gaussian: 1 point at the means; else Sobol').
  if (family == "gaussian") {
    n_int_eff <- 1L
    X_agd_int <- as.matrix(
      agd[, paste0(effect_modifiers, "_mean"), drop = FALSE])
    colnames(X_agd_int) <- NULL
  } else {
    n_int_eff <- as.integer(n_int)
    X_list <- lapply(seq_len(nrow(agd)), function(a) {
      means <- vapply(effect_modifiers,
                      function(v) agd[[paste0(v, "_mean")]][a], numeric(1))
      sds <- vapply(effect_modifiers,
                    function(v) agd[[paste0(v, "_sd")]][a], numeric(1))
      .cpaic_integration_points(means, sds, n_int_eff)
    })
    X_agd_int <- do.call(rbind, X_list)
  }

  base <- list(
    N_ipd = nrow(ipd), N_agd = nrow(agd), N_studies = length(studies),
    C = ncol(C), P = Q, Q = Q, n_int = n_int_eff,
    study_ipd = sidx(ipd[[study]]), Tc_ipd = Tc_ipd, X_ipd = X_ipd,
    em_idx = matrix(rep(seq_len(Q), each = nrow(ipd)), nrow = nrow(ipd)),
    study_agd = sidx(agd[[study]]), Tc_agd = Tc_agd, X_agd_int = X_agd_int,
    prior_intercept_sd = prior_intercept_sd, prior_beta_sd = prior_beta_sd,
    prior_reg_sd = prior_reg_sd
  )

  standata <- switch(
    family,
    binomial = c(base, list(
      y_ipd = as.integer(ipd[[outcome]]),
      r_agd = as.integer(agd[[r]]), n_agd = as.integer(agd[[n]]))),
    gaussian = c(base, list(
      y_ipd = as.numeric(ipd[[outcome]]),
      y_agd = as.numeric(agd[[outcome]]), se_agd = as.numeric(agd[[se]]))),
    poisson = c(base, list(
      y_ipd = as.integer(ipd[[outcome]]),
      offset_ipd = as.numeric(ipd[[exposure]]),
      r_agd = as.integer(agd[[r]]), E_agd = as.numeric(agd[[E]]))),
    survival = c(base, list(
      y_ipd = as.integer(ipd[[outcome]]),       # event indicator
      offset_ipd = as.numeric(ipd[[time]]),     # follow-up time
      r_agd = as.integer(agd[[r]]),             # total events
      E_agd = as.numeric(agd[[E]])))            # total person-time
  )

  mod <- .cpaic_stan_model(stan_family)
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
    row.names = NULL, stringsAsFactors = FALSE)

  structure(
    list(fit = fit, components = comp_tbl, C.matrix = C, comps = comps,
         family = family, effect_modifiers = effect_modifiers,
         method = "cML-NMR"),
    class = c("cpaic_mlnmr", "cpaic_fit"))
}

#' Compile (and cache) a cpaic Stan model with cmdstanr
#'
#' The model is compiled into a per-user cache directory rather than next to
#' the installed `.stan` file, so the package tree never acquires an
#' executable. The compiled binary is reused across calls.
#' @noRd
.cpaic_stan_model <- function(family) {
  stan_file <- system.file("stan", paste0("cpaic_", family, ".stan"),
                           package = "cpaic")
  if (stan_file == "") {
    stan_file <- file.path("inst", "stan", paste0("cpaic_", family, ".stan"))
  }
  cache <- tools::R_user_dir("cpaic", "cache")
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)
  dest <- file.path(cache, basename(stan_file))
  if (!file.exists(dest) || file.mtime(stan_file) > file.mtime(dest)) {
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
