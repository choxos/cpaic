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

#' Lexis (episode) expansion of survival IPD into interval-at-risk rows
#'
#' Splits each individual's follow-up at `cut_points`, returning one row per
#' interval the individual is at risk in, with the time at risk in that
#' interval (`.texp`), an event indicator for that interval (`.event`), and
#' the interval index (`.interval`). With `cut_points = NULL` there is a
#' single interval (the exponential model). Right-censoring only.
#' @noRd
.cpaic_lexis <- function(ipd, time_col, status_col, cut_points) {
  cuts <- c(0, cut_points, Inf)
  K <- length(cut_points) + 1L
  out <- vector("list", K)
  tt <- ipd[[time_col]]
  ss <- as.integer(ipd[[status_col]])
  for (k in seq_len(K)) {
    lo <- cuts[k]
    hi <- cuts[k + 1]
    at_risk <- tt > lo
    if (!any(at_risk)) next
    sub <- ipd[at_risk, , drop = FALSE]
    sub$.interval <- k
    sub$.texp <- pmin(tt[at_risk], hi) - lo
    sub$.event <- as.integer(ss[at_risk] == 1L & tt[at_risk] <= hi)
    out[[k]] <- sub
  }
  do.call(rbind, out)
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
#' @section Survival (approximations, read before use):
#' Survival uses a proportional-hazards model with a piecewise-constant
#' baseline log-hazard on the `cut_points` grid (`cut_points = NULL` gives the
#' exponential model). With `baseline = "mspline"` the interval heights are
#' *smoothed* by an M-spline basis evaluated at the interval midpoints. This is
#' a piecewise-exponential model with a smoothed step baseline; it is **not**
#' the continuous-time integrated M-spline survival likelihood of `multinma`,
#' which uses both the M-spline hazard basis and its integrated (I-spline)
#' cumulative hazard and supports left-, interval-censoring and delayed entry.
#' cpaic handles right-censoring only.
#'
#' Aggregate survival data are supplied as events and person-time per arm and
#' per interval, and the expected events are approximated by
#' `person-time x mean hazard`. This is an approximation even without
#' censoring, because higher-hazard individuals leave the risk set earlier, so
#' the *baseline* covariate distribution does not describe the covariate
#' distribution of the accumulated person-time. In a two-group example
#' (hazards 0.1 and 0.4, 50:50, follow-up to t = 10) it overstates expected
#' events by about 36%. Narrow intervals reduce the bias; supply IPD, or
#' interval-specific risk-set covariate summaries, where accuracy matters.
#'
#' @section Identifiability:
#' A relative effect is uniquely estimable only if its component contrast lies
#' in the row space of the within-study component design (Wigle et al. 2026);
#' [relative_effects()] returns `NA` otherwise rather than a prior-driven
#' number. Note this checks identification of `beta`; a component x
#' effect-modifier interaction is additionally identified only by covariate
#' variation on the contrasts that involve it, and interactions informed only
#' by aggregate arms are weakly identified (`prior_reg_sd` regularizes).
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
#' @param time,exposure IPD time / exposure column (survival, poisson).
#' @param r,n,E,se Aggregate columns: events `r`, sample size `n`
#'   (binomial), exposure `E` (poisson/survival), mean `outcome` and its
#'   standard error `se` (gaussian).
#' @param cut_points Survival only: interior interval boundaries for the
#'   piecewise-constant baseline. `NULL` (default) gives the exponential
#'   model; e.g. `c(6, 12)` gives three intervals.
#' @param interval Survival only: name of the aggregate interval-index
#'   column (values `1..K`) required when `cut_points` are supplied.
#' @param baseline Survival baseline hazard: `"piecewise"` (default, free
#'   step heights) or `"mspline"` (step heights smoothed by an M-spline
#'   evaluated at interval midpoints; see the Survival section).
#' @param n_basis Number of M-spline basis functions when
#'   `baseline = "mspline"` (must be `<=` the number of intervals).
#' @param cor Optional covariate correlation matrix for the Gaussian-copula
#'   integration. Must be a positive-definite correlation matrix (unit
#'   diagonal). Defaults to the within-study IPD correlation.
#' @param n_int Integration points per aggregate arm (ignored for
#'   `gaussian`, which is exact at the covariate means).
#' @param prior_intercept_sd,prior_beta_sd,prior_reg_sd Prior SDs.
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
                   r = "r", n = "n", E = "E", se = "se",
                   cut_points = NULL, interval = ".interval",
                   baseline = c("piecewise", "mspline"), n_basis = 6L,
                   cor = NULL, n_int = 64L,
                   prior_intercept_sd = 10, prior_beta_sd = 10,
                   prior_reg_sd = 2.5, chains = 4L, iter_warmup = 500L,
                   iter_sampling = 500L, seed = NULL, ...) {
  baseline <- match.arg(baseline)
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
    need(c(r, E), "agd")
    if (any(!pos(ipd[[time]]))) {
      stop("survival IPD time must be positive.", call. = FALSE)
    }
    if (!all(stats::na.omit(ipd[[outcome]]) %in% c(0, 1))) {
      stop("survival IPD status (`outcome`) must be 0 or 1.", call. = FALSE)
    }
    if (any(agd[[r]] < 0) || any(!pos(agd[[E]]))) {
      stop("aggregate survival needs `r` >= 0 and positive person-time `E`.",
           call. = FALSE)
    }
    warning("The aggregate survival likelihood approximates expected events ",
            "by person-time x mean hazard. Higher-hazard individuals leave ",
            "the risk set earlier, so this is biased upward unless the ",
            "intervals are narrow. See the Survival section of ?cmlnmr.",
            call. = FALSE)
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

  # First-order information design for the population-adjusted estimand
  # (beta, vec(Gamma)). This decides which contrasts are estimable AT A GIVEN
  # TARGET POPULATION, which is a strictly stronger requirement than
  # estimability of the component main effects. Built from the pre-expansion
  # IPD, since the Lexis expansion below duplicates rows.
  joint_design <- .cpaic_joint_design(C, ipd, agd, effect_modifiers,
                                      study = study, trt = trt)
  joint_design_ipd <- .cpaic_joint_design(C, ipd, agd[0, , drop = FALSE],
                                          effect_modifiers, study = study,
                                          trt = trt)

  # Survival: Lexis-expand IPD into interval-at-risk rows and configure the
  # piecewise baseline (K = 1 reduces to the exponential model).
  K <- 1L
  interval_ipd <- NULL
  time_ipd_vec <- NULL
  event_vec <- NULL
  Msp <- NULL
  interval_agd <- rep(1L, nrow(agd))
  if (family == "survival") {
    K <- length(cut_points) + 1L
    if (K > 1L) {
      if (!interval %in% names(agd)) {
        stop("With `cut_points`, `agd` must have an interval column `",
             interval, "` (values 1..", K, "), giving events `", r,
             "` and person-time `", E, "` per arm-interval.", call. = FALSE)
      }
      interval_agd <- as.integer(agd[[interval]])
      if (any(interval_agd < 1L | interval_agd > K)) {
        stop("`", interval, "` values must be in 1..", K, ".", call. = FALSE)
      }
    }
    if (baseline == "mspline") {
      if (K < 2L) {
        stop("baseline = 'mspline' needs `cut_points` to define a time grid.",
             call. = FALSE)
      }
      if (n_basis > K) {
        stop("`n_basis` (", n_basis, ") must be <= the number of intervals (",
             K, ").", call. = FALSE)
      }
      if (!requireNamespace("splines2", quietly = TRUE)) {
        stop("Package 'splines2' is required for baseline = 'mspline'.",
             call. = FALSE)
      }
      t_max <- max(c(ipd[[time]], cut_points))
      cuts_full <- c(0, cut_points, t_max)
      mid <- (utils::head(cuts_full, -1) + utils::tail(cuts_full, -1)) / 2
      Msp <- matrix(as.numeric(splines2::mSpline(
        mid, df = n_basis, degree = 3, intercept = TRUE,
        Boundary.knots = c(0, t_max))), nrow = K)
    }
    ipd <- .cpaic_lexis(ipd, time, outcome, cut_points)
    interval_ipd <- as.integer(ipd[[".interval"]])
    time_ipd_vec <- as.numeric(ipd[[".texp"]])
    event_vec <- as.integer(ipd[[".event"]])
  }

  Tc_ipd <- C[match(as.character(ipd[[trt]]), rownames(C)), , drop = FALSE]
  X_ipd <- as.matrix(ipd[, effect_modifiers, drop = FALSE])
  Tc_agd <- C[match(as.character(agd[[trt]]), rownames(C)), , drop = FALSE]

  # Integration points. Gaussian with all-normal margins is exact at the
  # covariate means (the identity link is linear in x), so a single point is
  # enough; any other case needs the QMC grid.
  gaussian_exact <- family == "gaussian" && all(margins == "normal")
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
      K = K,
      y_ipd = event_vec,                        # event in interval-row
      time_ipd = time_ipd_vec,                  # time at risk in interval
      interval_ipd = interval_ipd,
      r_agd = as.integer(agd[[r]]),             # events per arm-interval
      E_agd = as.numeric(agd[[E]]),             # person-time per arm-interval
      interval_agd = interval_agd))
  )
  if (family == "survival" && baseline == "mspline") {
    standata$n_basis <- ncol(Msp)
    standata$Msp <- Msp
  }

  mod <- .cpaic_stan_model(stan_family)
  fit <- mod$sample(
    data = standata, chains = chains, parallel_chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = seed %||% 1L, refresh = 0, show_messages = FALSE, ...)
  .cpaic_check_diagnostics(fit)

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

  structure(
    list(fit = fit, components = comp_tbl, C.matrix = C, comps = comps,
         family = family, effect_modifiers = effect_modifiers,
         margins = margins, cor = cor_mat, design = Dmat,
         null_space = null_space, rank = rk,
         joint_design = joint_design, joint_design_ipd = joint_design_ipd,
         reference = reference, sm = switch(family, binomial = "OR",
           gaussian = "MD", poisson = "IRR", survival = "HR"),
         method = "cML-NMR"),
    class = c("cpaic_mlnmr", "cpaic_fit"))
}

#' Warn on poor MCMC diagnostics from a cmdstanr fit
#' @noRd
.cpaic_check_diagnostics <- function(fit) {
  diag <- tryCatch(fit$diagnostic_summary(quiet = TRUE),
                   error = function(e) NULL)
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
  invisible(NULL)
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
component_effects.cpaic_mlnmr <- function(object, newdata = NULL, ...) {
  K <- ncol(object$C.matrix)
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
  # the component-main-effect design. Judging it on beta alone would report a
  # confident number for a component whose interaction is not identified.
  Id <- diag(K)
  V <- do.call(rbind, lapply(seq_len(K),
                             function(i) .cpaic_target_vec(Id[i, ], x)))
  bad <- !.cpaic_in_rowspace(V, .cpaic_null_space(object$joint_design))
  if (any(bad)) out[bad, -1L] <- NA_real_
  attr(out, "target") <- x
  out
}

#' @export
print.cpaic_mlnmr <- function(x, ...) {
  cat("cpaic: component-additive ML-NMR (Bayesian, ", x$family, ")\n",
      sep = "")
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
