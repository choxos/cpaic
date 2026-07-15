# Component-additive ML-NMR (Phase 2, Bayesian) ------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Build numerical integration points for an aggregate study
#'
#' Sobol' quasi-Monte-Carlo points are coupled by a Gaussian copula and then
#' pushed through each covariate's marginal inverse CDF, following the
#' construction used by `multinma::add_integration()`. Supported margins are
#' `"normal"` (from the study mean and SD) and `"bernoulli"` (from the study
#' mean, i.e. the prevalence). Using a normal margin for a binary covariate
#' generates integration points outside `{0, 1}` and integrates the
#' individual-level model over a population that does not exist, so the margin
#' must match the covariate.
#' @noRd
.cpaic_integration_points <- function(means, sds, n_int, cor = NULL,
                                      margins = NULL) {
  P <- length(means)
  if (!requireNamespace("randtoolbox", quietly = TRUE)) {
    stop("Package 'randtoolbox' is required for cmlnmr() integration.",
         call. = FALSE)
  }
  if (is.null(margins)) margins <- rep("normal", P)
  u <- matrix(randtoolbox::sobol(n = n_int, dim = P), nrow = n_int, ncol = P)
  # Guard the tails: qnorm(0)/qnorm(1) are infinite.
  eps <- 1 / (2 * max(n_int, 2L))
  u <- pmin(pmax(u, eps), 1 - eps)
  z <- stats::qnorm(u)                       # independent standard normals
  if (!is.null(cor) && P > 1L) z <- z %*% chol(cor)   # Gaussian copula

  out <- matrix(0, nrow = n_int, ncol = P)
  for (p in seq_len(P)) {
    out[, p] <- switch(
      margins[p],
      bernoulli = stats::qbinom(stats::pnorm(z[, p]), size = 1L,
                                prob = means[p]),
      normal = means[p] + sds[p] * z[, p],
      stop("Unsupported integration margin '", margins[p],
           "'; use \"normal\" or \"bernoulli\".", call. = FALSE)
    )
  }
  out
}

#' Guess a sensible integration margin for each effect modifier
#'
#' A covariate that is 0/1 in the IPD is treated as Bernoulli; anything else
#' as normal. Override with the `margins` argument of [cmlnmr()].
#' @noRd
.cpaic_guess_margins <- function(ipd, effect_modifiers) {
  vapply(effect_modifiers, function(v) {
    x <- stats::na.omit(ipd[[v]])
    if (length(x) && all(x %in% c(0, 1))) "bernoulli" else "normal"
  }, character(1))
}

#' Covariate correlation matrix for the Gaussian copula
#'
#' The copula needs the correlation *within* a population. Pooling all IPD
#' rows across studies and arms confounds the within-study association with
#' between-study shifts in the covariate means, which can manufacture a large
#' correlation where none exists. Correlations are therefore computed within
#' each IPD study and pooled on the Fisher z scale, weighted by `n - 3`, as in
#' `multinma`.
#' @noRd
.cpaic_copula_cor <- function(ipd, effect_modifiers, study_col, given = NULL) {
  Q <- length(effect_modifiers)
  if (Q < 2L) return(NULL)

  if (!is.null(given)) {
    R <- as.matrix(given)
    if (!is.matrix(R) || nrow(R) != Q || ncol(R) != Q || any(!is.finite(R)) ||
        !isTRUE(all.equal(R, t(R), check.attributes = FALSE))) {
      stop("`cor` must be a finite symmetric ", Q, "x", Q, " matrix.",
           call. = FALSE)
    }
    if (!isTRUE(all.equal(diag(R), rep(1, Q), check.attributes = FALSE))) {
      stop("`cor` must be a correlation matrix (unit diagonal).",
           call. = FALSE)
    }
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop("`cor` must be positive definite.", call. = FALSE)
    }
    return(R)
  }

  idx <- split(seq_len(nrow(ipd)), as.character(ipd[[study_col]]))
  Zsum <- matrix(0, Q, Q)
  Wsum <- matrix(0, Q, Q)
  for (rows in idx) {
    if (length(rows) < 4L) next
    Xs <- as.matrix(ipd[rows, effect_modifiers, drop = FALSE])
    R <- suppressWarnings(stats::cor(Xs, use = "pairwise.complete.obs"))
    R[!is.finite(R)] <- 0
    R <- pmin(pmax(R, -0.999), 0.999)
    w <- length(rows) - 3L
    Zsum <- Zsum + w * atanh(R)
    Wsum <- Wsum + w
  }
  if (all(Wsum == 0)) return(NULL)
  R <- tanh(Zsum / pmax(Wsum, 1))
  diag(R) <- 1
  R <- (R + t(R)) / 2
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  R
}

#' Component-additive multilevel network meta-regression (ML-NMR)
#'
#' The Bayesian flagship of cpaic. The relative effect of every treatment is
#' the sum of its component effects, estimated jointly from individual patient
#' data (IPD) and aggregate data (AgD). Aggregate arms are fitted by
#' integrating the individual-level model over each study's covariate
#' distribution, averaging the outcome on its natural scale (not the link
#' scale). Because disconnected sub-networks share component parameters, the
#' network is connected by construction.
#'
#' The model includes component x effect-modifier interactions `gamma`, so the
#' treatment effect is **population-specific**:
#' \deqn{\theta_t(x) = C_t' (\beta + \Gamma x).}
#' The component main effects `beta` are the effects at the covariate origin
#' (`x = 0`) and are *not* by themselves a population-adjusted quantity. Use
#' `newdata` in [relative_effects()] / [component_effects()] to obtain effects
#' in a named target population.
#'
#' Supported families: `"binomial"` (logit), `"gaussian"` (identity),
#' `"poisson"` (log), and `"survival"`.
#'
#' @section Integration:
#' Aggregate covariates are integrated with Sobol' quasi-Monte-Carlo points
#' coupled by a Gaussian copula, whose correlation is pooled *within* IPD
#' studies on the Fisher z scale (or supplied via `cor`). Each covariate is
#' pushed through its own marginal inverse CDF: `margins` may be `"normal"`
#' (using `x_mean` and `x_sd`) or `"bernoulli"` (using `x_mean` as the
#' prevalence). Margins default to Bernoulli for covariates that are 0/1 in the
#' IPD and normal otherwise; a normal margin on a binary covariate would
#' integrate over a population that cannot occur.
#'
#' @section Random effects:
#' `trt_effects = "random"` adds study-arm deviations around the
#' component-implied relative effects. Deviations use a non-centered
#' parameterization by default. Within a multi-arm study, deviations relative
#' to the study baseline have the standard NMA correlation of 0.5. The
#' heterogeneity standard deviation `tau` has a half-normal(0, 1) prior by
#' default. The centered parameterization is available only to reproduce
#' sampling comparisons.
#'
#' @section Priors:
#' Defaults follow the Stan prior-choice recommendations. Component effects
#' use normal(0, 2.5), component by effect-modifier interactions use
#' normal(0, 1), study intercepts use normal(0, 2.5), and `tau` uses
#' half-normal(0, 1). Interaction priors do real regularization when Gamma is
#' weakly identified, so every fitted object records the complete prior
#' specification. Use [prior_sensitivity()] to quantify contrast movement and
#' `prior_predictive = TRUE` with [prior_predictive_check()] to inspect prior
#' implications before fitting the likelihood.
#'
#' @section Survival:
#' Survival uses the exact individual likelihood ported from `multinma`
#' (Phillippo et al. 2020). The model evaluates a hazard basis and its
#' integrated cumulative-hazard basis at every outcome, interval start, and
#' delayed-entry time. It supports observed events, right censoring, left
#' censoring, interval censoring, and delayed entry. `baseline = "piecewise"`
#' gives a piecewise-exponential baseline; `baseline = "mspline"` gives a
#' continuous cubic M-spline baseline.
#'
#' Aggregate survival input must contain reconstructed event and censoring rows
#' with the same outcome-time columns as IPD, plus repeated arm-level covariate
#' summaries. The likelihood of every aggregate row is averaged over its
#' covariate integration points with `log_sum_exp`. Aggregate event counts and
#' person-time alone cannot recover this likelihood and are rejected explicitly.
#'
#' Two qualifications, so that "exact" is not read more broadly than it should
#' be.
#'
#' * **The likelihood is exact; the covariate integration is not.** Every
#'   individual contribution (event, right, left and interval censoring, delayed
#'   entry) is the exact analytic expression, verified against closed form to
#'   machine precision. The *aggregate* likelihood, however, averages that exact
#'   contribution over a finite quasi-Monte-Carlo grid of `n_int` covariate
#'   points, so it carries an integration error that shrinks with `n_int` but is
#'   not zero. Increase `n_int` and confirm that the estimates are stable before
#'   relying on them. The earlier person-time approximation was biased by 36% in
#'   a two-group example; that particular bias is removed, but a finite
#'   integration error remains.
#' * **Each study has its own baseline hazard shape.** Every study carries its
#'   own set of spline (or step) coefficients, smoothed toward a common shape by
#'   a shared random-walk scale (`prior_aux_sd`), so the treatment effects do not
#'   have to absorb baseline misfit. A single global spline basis is built from
#'   the pooled follow-up range, so a study with much shorter follow-up may not
#'   inform the coefficients of the latest basis functions; those are then
#'   determined by the smoothing prior rather than by that study's data.
#'
#' @section Identifiability:
#' A relative effect is uniquely estimable only if its component contrast lies
#' in the row space of the within-study component design (Wigle et al. 2026);
#' [relative_effects()] returns `NA` otherwise rather than a prior-driven
#' number. Note this checks identification of `beta`; a component x
#' effect-modifier interaction is additionally identified only by covariate
#' variation on the contrasts that involve it, and interactions informed only
#' by aggregate arms are weakly identified (`prior_gamma_scale` regularizes).
#'
#' @param ipd Individual patient data (one row per patient).
#' @param agd Aggregate data (one row per arm) with the per-study covariate
#'   summaries `x_mean` (and `x_sd` for normal margins) for each effect
#'   modifier `x`.
#' @param effect_modifiers Character vector of effect-modifier names.
#' @param margins Optional named character vector giving the integration
#'   margin of each effect modifier: `"normal"` or `"bernoulli"`. Defaults to
#'   Bernoulli for 0/1 covariates and normal otherwise.
#' @param inactive,sep.comps Component coding (see [cpaic_network()]).
#'   `inactive = NULL` gives the *unanchored* component parameterization, in
#'   which every unit receives its own parameter (Wigle & Béliveau 2022).
#' @param family One of `"binomial"`, `"gaussian"`, `"poisson"`,
#'   `"survival"`.
#' @param study,trt Column names (in both `ipd` and `agd`).
#' @param outcome IPD outcome column: 0/1 (binomial), numeric (gaussian),
#'   count (poisson), or the event indicator for survival.
#' @param time,exposure Outcome-time column for survival in both IPD and AgD;
#'   IPD exposure column for Poisson outcomes.
#' @param start,entry Survival columns giving the lower endpoint for
#'   interval-censored outcomes and the delayed-entry time. Missing columns
#'   imply zero.
#' @param r,n,E,se Aggregate columns: events `r`, sample size `n`
#'   (binomial), exposure `E` (poisson), mean `outcome` and its standard error
#'   `se` (gaussian).
#' @param cut_points Survival only: interior interval boundaries for a
#'   piecewise baseline. `NULL` gives the exponential model. This argument is
#'   ignored for a continuous M-spline baseline.
#' @param interval Retained for source compatibility; exact survival data do
#'   not use interval-indexed event counts.
#' @param baseline Survival baseline hazard: `"piecewise"` (default, free
#'   step heights) or `"mspline"` (a continuous cubic M-spline with its exact
#'   integrated basis).
#' @param n_basis Number of cubic M-spline basis functions. Must be at least 4.
#' @param cor Optional covariate correlation matrix for the Gaussian-copula
#'   integration. Must be a positive-definite correlation matrix (unit
#'   diagonal). Defaults to the within-study IPD correlation.
#' @param n_int Integration points per aggregate arm (ignored for `gaussian`,
#'   which is exact at the covariate means).
#'
#'   This is the main cost lever for the survival families. An aggregate survival
#'   arm is supplied as reconstructed pseudo-IPD, so the aggregate likelihood is
#'   evaluated once per (aggregate row x integration point): the work grows as
#'   `nrow(agd) * n_int`, and the default of 64 is expensive on a trial with
#'   several hundred reconstructed patients. Sampling is usually well behaved on
#'   the fixed-effects model (no divergences in the fixed-effects checks here),
#'   though the random-effects survival model can still produce a few divergent
#'   transitions and occasional rejected simplex proposals; inspect the
#'   diagnostics rather than assuming they are clean. If a survival
#'   fit is slow, reduce `n_int` before suspecting the geometry, and confirm the
#'   answer is stable with [plot_integration_error()].
#' @param QR Logical scalar. If `TRUE`, apply the scaled thin QR
#'   reparameterization used by `multinma` to the complete fixed-effects design
#'   matrix. This is only a reparameterization: it must not change the posterior
#'   distribution, only the geometry the sampler explores. The default is
#'   `FALSE`, matching `multinma`.
#'
#'   Do not turn this on expecting a free improvement. On the component networks
#'   tested here the fixed-effects design was not badly conditioned (a condition
#'   number near 19, in a network where every active treatment shared a
#'   component), and `QR = TRUE` gave *fewer* effective samples per second than
#'   `QR = FALSE`, with no divergent transitions either way. The intuition that a
#'   component design must be severely collinear, because one component recurs
#'   across many multi-component treatments, is not borne out: the study
#'   intercepts and the spread of the integration points keep the conditioning
#'   mild. Check `Z_cond` on the fit, and reach for `QR = TRUE` only when it is
#'   large.
#' @param prior_aux_sd Scale of the half-normal prior on the baseline-hazard
#'   smoothing parameter (survival families only). Each study has its own
#'   baseline hazard, given a first-order random-walk prior on the log spline
#'   coefficients with this shared smoothing scale. This is a simplified relative
#'   of the smoothing prior in `multinma`, not the same prior: smaller values
#'   shrink every study's baseline toward equal spline weights, which is a smooth
#'   default shape and not a constant hazard. The default of 1 follows the Stan
#'   recommendation of a half-normal(0, 1) prior for a hierarchical scale.
#' @param trt_effects Treatment-effect model: `"fixed"` or `"random"`.
#' @param re_parameterization Random-effects parameterization. The default
#'   `"noncentered"` should be used for inference; `"centered"` is provided for
#'   sampling diagnostics.
#' @param prior_intercept_sd,prior_beta_sd,prior_reg_sd Standard deviations for
#'   study-intercept, component-effect, and prognostic-regression normal priors.
#' @param prior_gamma_dist,prior_gamma_scale,prior_gamma_df Distribution, scale,
#'   and degrees of freedom for interaction priors. The Student t option uses
#'   the stated degrees of freedom.
#' @param prior_tau_dist,prior_tau_scale,prior_tau_df Distribution, scale, and
#'   degrees of freedom for the positive heterogeneity prior.
#' @param prior_predictive If `TRUE`, sample from the prior and omit the
#'   observed likelihood. Replicated outcomes remain available for
#'   [prior_predictive_check()].
#' @param chains,iter_warmup,iter_sampling,seed Passed to `cmdstanr`.
#' @param ... Passed to the `cmdstanr` sampler (e.g. `adapt_delta`).
#'
#' @return An object of class `cpaic_mlnmr` with the `cmdstanr` fit, the
#'   component design, and a tidy table of component effects.
#' @references
#' Phillippo DM, Dias S, Ades AE, et al. (2020). Multilevel network
#' meta-regression for population-adjusted treatment comparisons. *JRSS A*,
#' 183(3), 1189--1210.
#'
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @seealso [cmaic()], [cstc()], [cnma_bridge()], [estimable_effects()]
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
#' # Effects in a named target population (x1 = 0.2), not at the origin:
#' relative_effects(fit, newdata = data.frame(x1 = 0.2))
#' }
#' @export
cmlnmr <- function(ipd, agd, effect_modifiers, inactive = NULL,
                   sep.comps = "+", family = "binomial",
                   margins = NULL,
                   study = ".study", trt = ".trt", outcome = ".y",
                   time = ".time", exposure = ".exposure",
                   start = ".start", entry = ".entry",
                   r = "r", n = "n", E = "E", se = "se",
                   cut_points = NULL, interval = ".interval",
                   baseline = c("piecewise", "mspline"), n_basis = 6L,
                   cor = NULL, n_int = 64L, QR = FALSE,
                   trt_effects = c("fixed", "random"),
                   re_parameterization = c("noncentered", "centered"),
                   prior_intercept_sd = 2.5, prior_aux_sd = 1,
                   prior_beta_sd = 2.5,
                   prior_reg_sd = 1,
                   prior_gamma_dist = c("normal", "student_t"),
                   prior_gamma_scale = 1, prior_gamma_df = 4,
                   prior_tau_dist = c("half-normal", "half-student-t"),
                   prior_tau_scale = 1, prior_tau_df = 4,
                   prior_predictive = FALSE,
                   chains = 4L, iter_warmup = 500L,
                   iter_sampling = 500L, seed = NULL, ...) {
  baseline <- match.arg(baseline)
  trt_effects <- match.arg(trt_effects)
  re_parameterization <- match.arg(re_parameterization)
  prior_gamma_dist <- match.arg(prior_gamma_dist)
  prior_tau_dist <- match.arg(prior_tau_dist)
  family <- match.arg(family,
                      c("binomial", "gaussian", "poisson", "survival"))
  stan_family <- switch(family, binomial = "binomial", gaussian = "normal",
                        poisson = "poisson", survival = "survival")
  if (family == "survival" && baseline == "mspline") {
    stan_family <- "survival_mspline"
  }
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmlnmr() needs 'cmdstanr'. Install from ",
         "https://stan-dev.r-universe.dev .", call. = FALSE)
  }
  if (!is.logical(prior_predictive) || length(prior_predictive) != 1L ||
      is.na(prior_predictive)) {
    stop("`prior_predictive` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(QR) || length(QR) != 1L || is.na(QR)) {
    stop("`QR` must be TRUE or FALSE.", call. = FALSE)
  }
  prior_values <- c(
    prior_intercept_sd = prior_intercept_sd,
    prior_aux_sd = prior_aux_sd,
    prior_beta_sd = prior_beta_sd,
    prior_reg_sd = prior_reg_sd,
    prior_gamma_scale = prior_gamma_scale,
    prior_gamma_df = prior_gamma_df,
    prior_tau_scale = prior_tau_scale,
    prior_tau_df = prior_tau_df
  )
  if (any(!is.finite(prior_values)) || any(prior_values <= 0)) {
    stop("All prior scales and degrees of freedom must be positive and finite.",
         call. = FALSE)
  }

  sample_args <- list(...)
  if (length(sample_args) &&
      (is.null(names(sample_args)) || any(!nzchar(names(sample_args))))) {
    stop("Sampler arguments in `...` must be named.", call. = FALSE)
  }
  ipd_input <- as.data.frame(ipd)
  agd_input <- as.data.frame(agd)
  ipd <- ipd_input
  agd <- agd_input
  Q <- length(effect_modifiers)

  all_trts <- sort(unique(c(as.character(ipd[[trt]]),
                            as.character(agd[[trt]]))))
  C <- build_C_matrix(all_trts, sep.comps = sep.comps, inactive = inactive)
  comps <- colnames(C)
  studies <- sort(unique(c(as.character(ipd[[study]]),
                           as.character(agd[[study]]))))
  sidx <- function(s) match(as.character(s), studies)

  # Front-door validation.
  if (Q < 1L) {
    stop("`effect_modifiers` must name at least one covariate.", call. = FALSE)
  }
  if (nrow(ipd) == 0L) {
    stop("cmlnmr() requires at least one row of IPD in this version.",
         call. = FALSE)
  }
  if (nrow(agd) == 0L) {
    stop("cmlnmr() requires aggregate data.", call. = FALSE)
  }
  if (!is.null(inactive) && !inactive %in% all_trts) {
    stop("`inactive` (\"", inactive, "\") is not one of the treatments.",
         call. = FALSE)
  }

  # Every study must be at least two-arm. A single-arm study carries no contrast,
  # so it can contribute nothing to a relative effect through the linear
  # predictor; its study intercept absorbs its data entirely. Whatever it appears
  # to contribute comes only from the curvature of the aggregate integral, which
  # is precisely the weak, ecological channel that should not be relied on
  # silently. Reject it rather than let it influence the treatment effects.
  arms_per_study <- tapply(
    c(as.character(ipd[[trt]]), as.character(agd[[trt]])),
    c(as.character(ipd[[study]]), as.character(agd[[study]])),
    function(z) length(unique(z)))
  lone <- names(arms_per_study)[arms_per_study < 2L]
  if (length(lone)) {
    stop("Single-arm study/studies: ", paste(lone, collapse = ", "),
         ". A study with one arm carries no treatment contrast and cannot inform ",
         "a relative effect; its study intercept absorbs it entirely. Remove it, ",
         "or supply its comparator arm.", call. = FALSE)
  }
  if (length(cut_points) &&
      (is.unsorted(cut_points, strictly = TRUE) || any(cut_points <= 0))) {
    stop("`cut_points` must be strictly increasing and positive.",
         call. = FALSE)
  }
  if (as.integer(n_int) < 1L) {
    stop("`n_int` must be a positive integer.", call. = FALSE)
  }
  miss_em <- setdiff(effect_modifiers, names(ipd))
  if (length(miss_em)) {
    stop("`ipd` is missing effect-modifier column(s): ",
         paste(miss_em, collapse = ", "), call. = FALSE)
  }

  # Integration margins: Bernoulli for 0/1 covariates unless overridden.
  guessed <- .cpaic_guess_margins(ipd, effect_modifiers)
  if (is.null(margins)) {
    margins <- guessed
  } else {
    if (is.null(names(margins))) names(margins) <- effect_modifiers
    unknown <- setdiff(names(margins), effect_modifiers)
    if (length(unknown)) {
      stop("`margins` names unknown effect modifier(s): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    bad <- setdiff(margins, c("normal", "bernoulli"))
    if (length(bad)) {
      stop("`margins` must be \"normal\" or \"bernoulli\"; got: ",
           paste(unique(bad), collapse = ", "), call. = FALSE)
    }
    guessed[names(margins)] <- margins
    margins <- guessed
  }
  margins <- margins[effect_modifiers]

  for (v in effect_modifiers) {
    if (!paste0(v, "_mean") %in% names(agd)) {
      stop("`agd` is missing `", v, "_mean`.", call. = FALSE)
    }
    if (margins[[v]] == "bernoulli") {
      p <- agd[[paste0(v, "_mean")]]
      if (any(!is.finite(p) | p < 0 | p > 1)) {
        stop("`", v, "_mean` is a Bernoulli prevalence and must lie in [0, 1].",
             call. = FALSE)
      }
    } else if (family != "gaussian") {
      sdn <- paste0(v, "_sd")
      if (!sdn %in% names(agd)) stop("`agd` is missing `", sdn, "`.",
                                     call. = FALSE)
      if (any(!is.finite(agd[[sdn]]) | agd[[sdn]] <= 0)) {
        stop("`", sdn, "` must be positive and finite.", call. = FALSE)
      }
    }
  }

  # Family-specific value checks (guard against invalid Stan inputs).
  surv_ipd <- NULL
  surv_agd <- NULL
  need <- function(cols, where) {
    m <- setdiff(cols, names(if (where == "agd") agd else ipd))
    if (length(m)) stop("`", where, "` is missing column(s): ",
                        paste(m, collapse = ", "), call. = FALSE)
  }
  pos <- function(x) is.finite(x) & x > 0
  if (family == "binomial") {
    need(c(r, n), "agd")
    if (!all(stats::na.omit(ipd[[outcome]]) %in% c(0, 1))) {
      stop("binomial IPD outcome must be 0 or 1.", call. = FALSE)
    }
    if (any(agd[[r]] < 0 | agd[[r]] > agd[[n]])) {
      stop("aggregate binomial counts must satisfy 0 <= r <= n.",
           call. = FALSE)
    }
  } else if (family == "gaussian") {
    need(c(outcome, se), "agd")
    if (any(!pos(agd[[se]]))) {
      stop("aggregate gaussian `se` must be positive.", call. = FALSE)
    }
  } else if (family == "poisson") {
    need(c(r, E), "agd")
    if (any(agd[[r]] < 0) || any(!pos(agd[[E]]))) {
      stop("aggregate poisson needs `r` >= 0 and positive exposure `E`.",
           call. = FALSE)
    }
    if (!is.null(exposure) && exposure %in% names(ipd) &&
        any(!pos(ipd[[exposure]]))) {
      stop("poisson IPD exposure must be positive.", call. = FALSE)
    }
  } else if (family == "survival") {
    surv_ipd <- .cpaic_survival_outcomes(
      ipd, time_col = time, status_col = outcome, start_col = start,
      entry_col = entry, source = "ipd"
    )
    surv_agd <- .cpaic_survival_outcomes(
      agd, time_col = time, status_col = outcome, start_col = start,
      entry_col = entry, source = "agd"
    )
  }

  # Reference treatment for reported relative effects.
  reference <- if (!is.null(inactive) && inactive %in% rownames(C)) {
    inactive
  } else {
    rownames(C)[1]
  }

  # Estimability: the within-study arm contrasts, in component space, form the
  # design whose row space contains exactly the estimable relative effects.
  arm_tbl <- unique(rbind(
    data.frame(s = as.character(ipd[[study]]), t = as.character(ipd[[trt]]),
               stringsAsFactors = FALSE),
    data.frame(s = as.character(agd[[study]]), t = as.character(agd[[trt]]),
               stringsAsFactors = FALSE)))
  Dlist <- lapply(unique(arm_tbl$s), function(ss) {
    ts <- arm_tbl$t[arm_tbl$s == ss]
    if (length(ts) < 2L) return(NULL)
    Cs <- C[match(ts, rownames(C)), , drop = FALSE]
    sweep(Cs[-1, , drop = FALSE], 2, Cs[1, ], "-")
  })
  Dmat <- do.call(rbind, Dlist)
  if (is.null(Dmat) || !nrow(Dmat)) {
    Dmat <- matrix(0, nrow = 1L, ncol = ncol(C))
  }
  rk <- as.integer(Matrix::rankMatrix(Dmat))
  null_space <- .cpaic_null_space(Dmat)
  if (rk < ncol(C)) {
    est_tbl <- .cpaic_estimable_table(C, null_space, reference)
    bad <- est_tbl$treatment[!est_tbl$estimable]
    warning("Only ", rk, " of ", ncol(C), " component effects are identified ",
            "by within-study contrasts; the rest are informed only by the ",
            "prior. ",
            if (length(bad))
              paste0("Relative effects versus \"", reference,
                     "\" that are NOT estimable (returned as NA): ",
                     paste(bad, collapse = ", "), ".")
            else "All relative effects versus the reference remain estimable.",
            call. = FALSE)
  }

  # Covariate correlation for the Gaussian-copula integration, pooled WITHIN
  # IPD studies (pooling across studies would confound within-study
  # association with between-study mean shifts).
  cor_mat <- .cpaic_copula_cor(ipd, effect_modifiers, study, cor)

  # The copula takes its correlation on the LATENT Gaussian scale. When the
  # correlation is estimated automatically it is an observed Pearson correlation,
  # which equals the latent one only for normal margins. For a Bernoulli margin
  # the two differ, so an auto-estimated correlation reproduces the observed
  # association only approximately. A correct observed-to-latent transform (as in
  # multinma) is not applied here; supply `cor` on the latent scale if the exact
  # association matters. Only warn when it can actually bite: two or more
  # correlated modifiers with at least one Bernoulli margin, and no user `cor`.
  if (is.null(cor) && !is.null(cor_mat) && any(margins == "bernoulli") &&
      length(effect_modifiers) >= 2L) {
    warning("The auto-estimated covariate correlation is an observed Pearson ",
            "correlation but is used as a latent copula correlation. For the ",
            "Bernoulli margin(s) present these differ, so the integration ",
            "reproduces the observed association only approximately. Supply ",
            "`cor` on the latent scale for an exact match.", call. = FALSE)
  }

  # First-order information design for the population-adjusted estimand
  # (beta, vec(Gamma)). This decides which contrasts are estimable AT A GIVEN
  # TARGET POPULATION, which is a strictly stronger requirement than
  # estimability of the component main effects.
  joint_design <- .cpaic_joint_design(C, ipd, agd, effect_modifiers,
                                      study = study, trt = trt)
  joint_design_ipd <- .cpaic_joint_design(C, ipd, agd[0, , drop = FALSE],
                                          effect_modifiers, study = study,
                                          trt = trt)

  # Exact survival bases. The integrated basis is evaluated at every outcome,
  # interval-start, and delayed-entry time.
  survival_spec <- NULL
  survival_ipd_basis <- NULL
  survival_agd_basis <- NULL
  if (family == "survival") {
    basis_cuts <- cut_points
    if (baseline == "mspline" && length(cut_points)) {
      warning("`cut_points` is ignored for baseline = 'mspline'; `n_basis` ",
              "controls the continuous cubic basis.", call. = FALSE)
      basis_cuts <- numeric()
    }
    survival_spec <- .cpaic_survival_basis_spec(
      observed_times = c(surv_ipd$time, surv_agd$time),
      baseline = baseline, cut_points = basis_cuts, n_basis = n_basis
    )
    survival_ipd_basis <- .cpaic_survival_basis_eval(
      survival_spec, surv_ipd$time, surv_ipd$start_time,
      surv_ipd$entry_time
    )
    survival_agd_basis <- .cpaic_survival_basis_eval(
      survival_spec, surv_agd$time, surv_agd$start_time,
      surv_agd$entry_time
    )
  }

  Tc_ipd <- C[match(as.character(ipd[[trt]]), rownames(C)), , drop = FALSE]
  X_ipd <- as.matrix(ipd[, effect_modifiers, drop = FALSE])
  Tc_agd <- C[match(as.character(agd[[trt]]), rownames(C)), , drop = FALSE]

  # Integration points. Under the identity link the aggregate mean is a linear
  # functional of the covariates, so it depends only on their means: E[eta] =
  # eta(E[x]), whatever the marginal distributions or their correlation. A single
  # point at the covariate means is therefore EXACT for the gaussian family
  # regardless of margin type, including Bernoulli. (This path already evaluates
  # X_agd_int at the *_mean columns, which is E[x] for a Bernoulli prevalence as
  # well as for a normal mean.) Restricting the exact path to all-normal margins
  # would force a rare binary modifier through finite QMC and shift its mean; the
  # nonlinear families still need the QMC grid.
  gaussian_exact <- family == "gaussian"
  if (gaussian_exact) {
    n_int_eff <- 1L
    X_agd_int <- as.matrix(
      agd[, paste0(effect_modifiers, "_mean"), drop = FALSE])
    colnames(X_agd_int) <- NULL
  } else {
    n_int_eff <- as.integer(n_int)
    sd_of <- function(v, a) {
      if (margins[[v]] == "bernoulli") return(NA_real_)
      agd[[paste0(v, "_sd")]][a]
    }
    X_list <- lapply(seq_len(nrow(agd)), function(a) {
      means <- vapply(effect_modifiers,
                      function(v) agd[[paste0(v, "_mean")]][a], numeric(1))
      sds <- vapply(effect_modifiers, sd_of, numeric(1), a = a)
      .cpaic_integration_points(means, sds, n_int_eff, cor = cor_mat,
                                margins = margins)
    })
    X_agd_int <- do.call(rbind, X_list)
  }

  re <- .cpaic_random_effects(ipd, agd, study, trt)
  if (trt_effects == "random" && re$N_delta < 1L) {
    stop("Random effects need at least one within-study treatment contrast.",
         call. = FALSE)
  }
  N_delta_stan <- max(1L, re$N_delta)
  L_delta_stan <- if (re$N_delta) re$L_delta else matrix(1, 1, 1)

  study_ipd <- sidx(ipd[[study]])
  study_agd <- sidx(agd[[study]])
  agd_int_idx <- rep(seq_len(nrow(agd)), each = n_int_eff)
  Z_ipd <- .cpaic_fixed_design(
    study = study_ipd, Tc = Tc_ipd, X = X_ipd,
    n_studies = length(studies), emc = seq_len(Q)
  )
  Z_agd_int <- .cpaic_fixed_design(
    study = rep(study_agd, each = n_int_eff),
    Tc = Tc_agd[agd_int_idx, , drop = FALSE], X = X_agd_int,
    n_studies = length(studies), emc = seq_len(Q)
  )
  Z <- rbind(Z_ipd, Z_agd_int)
  # A component design is collinear by construction: a component that appears in
  # several multi-component treatments has a column correlated with each of
  # them. The condition number says how badly, and so whether `QR = TRUE` is
  # likely to help the sampler.
  Z_cond <- tryCatch(kappa(Z), error = function(e) NA_real_)
  if (QR) {
    qr_design <- .cpaic_thin_qr(Z)
    Z_stan <- qr_design$Q
    R_inv <- qr_design$R_inv
  } else {
    Z_stan <- Z
    R_inv <- matrix(0, nrow = 0L, ncol = 0L)
  }
  Z_ipd_stan <- Z_stan[seq_len(nrow(ipd)), , drop = FALSE]
  Z_agd_int_stan <- Z_stan[nrow(ipd) + seq_len(nrow(X_agd_int)), ,
                            drop = FALSE]

  base <- list(
    N_ipd = nrow(ipd), N_agd = nrow(agd), N_studies = length(studies),
    C = ncol(C), P = Q, Q = Q, n_int = n_int_eff,
    nX = ncol(Z), QR = as.integer(QR),
    Z_ipd = Z_ipd_stan, Z_agd_int = Z_agd_int_stan, R_inv = R_inv,
    RE = as.integer(trt_effects == "random"),
    noncentered = as.integer(re_parameterization == "noncentered"),
    N_delta = N_delta_stan, L_delta = L_delta_stan,
    re_idx_ipd = re$re_idx_ipd, re_idx_agd = re$re_idx_agd,
    prior_only = as.integer(prior_predictive),
    prior_intercept_sd = prior_intercept_sd, prior_beta_sd = prior_beta_sd,
    prior_aux_sd = prior_aux_sd, prior_reg_sd = prior_reg_sd,
    prior_gamma_dist = match(prior_gamma_dist, c("normal", "student_t")),
    prior_gamma_scale = prior_gamma_scale, prior_gamma_df = prior_gamma_df,
    prior_tau_dist = match(prior_tau_dist,
                           c("half-normal", "half-student-t")),
    prior_tau_scale = prior_tau_scale, prior_tau_df = prior_tau_df
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
      N_base = survival_spec$n_basis,
      # The study index is carried separately for the survival families. The
      # study intercept is inside the QR design, but each study also has its OWN
      # baseline hazard, which has to be indexed by study.
      study_ipd = study_ipd, study_agd = study_agd,
      time_basis_ipd = survival_ipd_basis$time,
      itime_basis_ipd = survival_ipd_basis$itime,
      start_basis_ipd = survival_ipd_basis$start_itime,
      entry_basis_ipd = survival_ipd_basis$entry_itime,
      delayed_ipd = survival_ipd_basis$delayed,
      status_ipd = surv_ipd$status,
      time_basis_agd = survival_agd_basis$time,
      itime_basis_agd = survival_agd_basis$itime,
      start_basis_agd = survival_agd_basis$start_itime,
      entry_basis_agd = survival_agd_basis$entry_itime,
      delayed_agd = survival_agd_basis$delayed,
      status_agd = surv_agd$status))
  )

  mod <- .cpaic_stan_model(stan_family)
  sample_defaults <- list(
    data = standata, chains = chains, parallel_chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed %||% 1L, refresh = 0, show_messages = FALSE
  )
  if (length(sample_args)) {
    sample_defaults[names(sample_args)] <- NULL
  }
  fit <- do.call(mod$sample, c(sample_defaults, sample_args))
  diagnostics <- .cpaic_check_diagnostics(fit)

  beta_draws <- fit$draws("beta", format = "draws_matrix")
  comp_tbl <- data.frame(
    component = comps,
    estimate = apply(beta_draws, 2, mean),
    se = apply(beta_draws, 2, stats::sd),
    lower = apply(beta_draws, 2, stats::quantile, 0.025),
    upper = apply(beta_draws, 2, stats::quantile, 0.975),
    row.names = NULL, stringsAsFactors = FALSE)
  bad_comp <- !.cpaic_in_rowspace(diag(ncol(C)), null_space)
  if (any(bad_comp)) {
    comp_tbl[bad_comp, c("estimate", "se", "lower", "upper")] <- NA_real_
  }

  observed <- switch(
    family,
    binomial = list(ipd = as.integer(ipd[[outcome]]),
                    agd = as.integer(agd[[r]])),
    gaussian = list(ipd = as.numeric(ipd[[outcome]]),
                    agd = as.numeric(agd[[outcome]])),
    poisson = list(ipd = as.integer(ipd[[outcome]]),
                   agd = as.integer(agd[[r]])),
    survival = list(ipd = as.integer(surv_ipd$status != 0L),
                    agd = as.integer(surv_agd$status != 0L))
  )
  rep_variables <- switch(
    family,
    binomial = list(ipd = "yrep_ipd", agd = "rrep_agd"),
    gaussian = list(ipd = "yrep_ipd", agd = "yrep_agd"),
    poisson = list(ipd = "yrep_ipd", agd = "rrep_agd"),
    survival = list(ipd = "event_rep_ipd", agd = "event_rep_agd")
  )
  refit_args <- list(
    ipd = ipd_input, agd = agd_input,
    effect_modifiers = effect_modifiers, inactive = inactive,
    sep.comps = sep.comps, family = family, margins = margins,
    study = study, trt = trt, outcome = outcome, time = time,
    exposure = exposure, start = start, entry = entry,
    r = r, n = n, E = E, se = se, cut_points = cut_points,
    interval = interval, baseline = baseline, n_basis = n_basis,
    cor = cor, n_int = n_int, QR = QR, trt_effects = trt_effects,
    re_parameterization = re_parameterization,
    prior_intercept_sd = prior_intercept_sd,
    prior_aux_sd = prior_aux_sd,
    prior_beta_sd = prior_beta_sd, prior_reg_sd = prior_reg_sd,
    prior_gamma_dist = prior_gamma_dist,
    prior_gamma_scale = prior_gamma_scale, prior_gamma_df = prior_gamma_df,
    prior_tau_dist = prior_tau_dist, prior_tau_scale = prior_tau_scale,
    prior_tau_df = prior_tau_df, prior_predictive = prior_predictive,
    chains = chains, iter_warmup = iter_warmup,
    iter_sampling = iter_sampling, seed = seed
  )
  if (length(sample_args)) refit_args[names(sample_args)] <- sample_args

  structure(
    list(fit = fit, components = comp_tbl, C.matrix = C, comps = comps,
         family = family, effect_modifiers = effect_modifiers,
         margins = margins, cor = cor_mat, QR = QR, Z_cond = Z_cond,
         design = Dmat,
         null_space = null_space, rank = rk,
         joint_design = joint_design, joint_design_ipd = joint_design_ipd,
         reference = reference, sm = switch(family, binomial = "OR",
           gaussian = "MD", poisson = "IRR", survival = "HR"),
         method = "cML-NMR", trt_effects = trt_effects,
         re_parameterization = re_parameterization,
         priors = list(
           intercept = list(distribution = "normal", location = 0,
                            scale = prior_intercept_sd),
           beta = list(distribution = "normal", location = 0,
                       scale = prior_beta_sd),
           regression = list(distribution = "normal", location = 0,
                             scale = prior_reg_sd),
           gamma = list(distribution = prior_gamma_dist, location = 0,
                        scale = prior_gamma_scale, df = prior_gamma_df),
           tau = list(distribution = prior_tau_dist, location = 0,
                      scale = prior_tau_scale, df = prior_tau_df)),
         prior_predictive = prior_predictive, observed = observed,
         rep_variables = rep_variables, diagnostics = diagnostics,
         refit_args = refit_args),
    class = c("cpaic_mlnmr", "cpaic_fit"))
}

#' Warn on poor MCMC diagnostics from a cmdstanr fit
#' @noRd
.cpaic_check_diagnostics <- function(fit) {
  diag <- tryCatch(fit$diagnostic_summary(quiet = TRUE),
                   error = function(e) NULL)
  nd <- 0L
  ntd <- 0L
  if (!is.null(diag)) {
    nd <- sum(diag$num_divergent %||% 0)
    if (nd > 0) {
      warning(nd, " divergent transition(s) in cmlnmr(); results may be ",
              "unreliable (consider higher adapt_delta or more iterations).",
              call. = FALSE)
    }
    ntd <- sum(diag$num_max_treedepth %||% 0)
    if (ntd > 0) {
      warning(ntd, " iteration(s) saturated the maximum tree depth.",
              call. = FALSE)
    }
  }
  rh <- tryCatch(max(fit$summary(c("beta", "mu"))$rhat, na.rm = TRUE),
                 error = function(e) NA_real_)
  if (is.finite(rh) && rh > 1.05) {
    warning("Maximum Rhat = ", round(rh, 3), " (> 1.05); the chains may not ",
            "have converged. Increase `iter_warmup`/`iter_sampling`.",
            call. = FALSE)
  }
  invisible(list(divergences = as.integer(nd),
                 max_treedepth = as.integer(ntd), max_rhat = rh))
}

#' Compile (and cache) a cpaic Stan model with cmdstanr
#'
#' The model is compiled into a per-user cache directory rather than next to the
#' installed `.stan` file, so the package tree never acquires an executable. The
#' compiled binary is reused across calls.
#'
#' The cache is keyed by the CONTENT of the Stan source, not by its file name.
#' Keying on the name and refreshing on the modification time is not safe:
#' installing a package does not reliably give the installed `.stan` file a
#' modification time later than a binary already sitting in the cache, so an
#' upgraded cpaic could silently keep running the previous version's compiled
#' model. A content hash cannot go stale, and it also lets two versions of the
#' package coexist, which is what makes the QR equivalence test possible.
#' @noRd
.cpaic_stan_model <- function(family) {
  stan_file <- system.file("stan", paste0("cpaic_", family, ".stan"),
                           package = "cpaic")
  if (stan_file == "") {
    stan_file <- file.path("inst", "stan", paste0("cpaic_", family, ".stan"))
  }
  cache <- tools::R_user_dir("cpaic", "cache")
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE)

  key <- unname(tools::md5sum(stan_file))
  if (is.na(key)) {
    stop("Could not find the Stan model for family '", family, "'.",
         call. = FALSE)
  }
  dest <- file.path(cache, sprintf("cpaic_%s_%s.stan", family,
                                   substr(key, 1L, 10L)))
  if (!file.exists(dest)) file.copy(stan_file, dest, overwrite = TRUE)
  cmdstanr::cmdstan_model(dest)
}

#' @export
component_effects.cpaic_mlnmr <- function(object, newdata = NULL, ...) {
  Q <- length(object$effect_modifiers)
  x <- if (is.null(newdata)) rep(0, Q) else
    .cpaic_target_x(newdata, object$effect_modifiers)
  Beff <- .cpaic_beta_at(object, x)
  out <- data.frame(
    component = object$comps,
    estimate = apply(Beff, 2, mean),
    se = apply(Beff, 2, stats::sd),
    lower = apply(Beff, 2, stats::quantile, 0.025),
    upper = apply(Beff, 2, stats::quantile, 0.975),
    row.names = NULL, stringsAsFactors = FALSE)
  # A component effect IN A TARGET POPULATION is beta_c + Gamma_c' x, so its
  # estimability must be judged against the JOINT design at that x, not against
  # the component-main-effect design alone.
  Id <- diag(ncol(object$C.matrix))
  Vv <- do.call(rbind, lapply(seq_len(nrow(Id)),
                              function(i) .cpaic_target_vec(Id[i, ], x)))
  bad <- !.cpaic_in_rowspace(Vv, .cpaic_null_space(object$joint_design))
  if (any(bad)) out[bad, -1L] <- NA_real_
  attr(out, "target") <- x
  out
}

#' @export
print.cpaic_mlnmr <- function(x, ...) {
  cat("cpaic: component-additive ML-NMR (Bayesian, ", x$family, ")\n",
      sep = "")
  cat("  Treatment effects: ", x$trt_effects,
      if (x$trt_effects == "random") {
        paste0(" (", x$re_parameterization, ")")
      } else "", "\n", sep = "")
  if (isTRUE(x$prior_predictive)) {
    cat("  Prior-predictive fit: observed likelihood omitted.\n")
  }
  cat("  Effect modifiers: ",
      paste(names(x$margins), " [", x$margins, "]", sep = "",
            collapse = ", "), "\n", sep = "")
  cat("  Component effects below are at the covariate origin (x = 0).\n")
  cat("  For a target population use relative_effects(fit, newdata = ...).\n\n")
  comp <- x$components
  comp[, c("estimate", "se", "lower", "upper")] <-
    round(comp[, c("estimate", "se", "lower", "upper")], 3)
  print(comp, row.names = FALSE)
  invisible(x)
}
