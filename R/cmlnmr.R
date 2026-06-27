# Component-additive ML-NMR (Phase 2, Bayesian) ------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Build numerical integration points for an aggregate study
#'
#' Sobol' quasi-Monte-Carlo points mapped through normal margins defined by
#' the study covariate means and standard deviations. When a correlation
#' matrix `cor` is supplied the margins are coupled by a Gaussian copula
#' (for normal margins this is exactly correlated multivariate normal draws,
#' `mean + sd * z` with `z` correlated by `chol(cor)`); otherwise the
#' covariates are integrated independently.
#' @noRd
.cpaic_integration_points <- function(means, sds, n_int, cor = NULL) {
  P <- length(means)
  if (!requireNamespace("randtoolbox", quietly = TRUE)) {
    stop("Package 'randtoolbox' is required for cmlnmr() integration.",
         call. = FALSE)
  }
  u <- matrix(randtoolbox::sobol(n = n_int, dim = P), nrow = n_int, ncol = P)
  z <- stats::qnorm(u)                       # independent standard normals
  if (!is.null(cor) && P > 1L) z <- z %*% chol(cor)   # Gaussian copula
  sweep(sweep(z, 2, sds, "*"), 2, means, "+")
}

#' Covariate correlation matrix for the Gaussian copula, from IPD
#' @noRd
.cpaic_copula_cor <- function(X, given = NULL) {
  if (ncol(X) < 2L) return(NULL)
  R <- if (is.null(given)) stats::cor(X, use = "pairwise.complete.obs") else given
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
#' single interval (the exponential model).
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
#' the sum of its component effects (`theta = C beta`), estimated jointly
#' from individual patient data (IPD) and aggregate data (AgD). Aggregate
#' arms are fitted by integrating the individual-level model over each
#' study's covariate distribution. Because disconnected sub-networks share
#' component parameters, the network is connected by construction.
#'
#' Supported families: `"binomial"` (logit), `"gaussian"` (identity),
#' `"poisson"` (log), and `"survival"`. Survival uses a proportional-hazards
#' model with a flexible baseline: `baseline = "piecewise"` gives a
#' piecewise-exponential step baseline (one level per interval defined by
#' `cut_points`; exponential when `cut_points = NULL`), while
#' `baseline = "mspline"` gives a smooth baseline hazard from a non-negative
#' M-spline basis (`n_basis` functions, simplex weights) evaluated at the
#' interval midpoints over the `cut_points` grid. Individual patient data are
#' split at the cut points internally; aggregate survival data are supplied
#' as events and person-time per arm and (when `cut_points` are given) per
#' interval. The aggregate likelihood approximates the expected events in an
#' arm-interval by person-time times the integrated hazard; this assumes the
#' person-time is independent of the covariates within an interval (most
#' accurate when the intervals are narrow), so it is approximate when effect
#' modifiers also drive censoring.
#'
#' Aggregate covariates are integrated with a Gaussian copula: the
#' correlation is estimated from the IPD (or supplied via `cor`) so the
#' integration points respect the covariate correlation structure rather
#' than treating the effect modifiers as independent.
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
#' @param cut_points Survival only: interior interval boundaries for the
#'   piecewise-exponential baseline. `NULL` (default) gives the exponential
#'   model; e.g. `c(6, 12)` gives three intervals.
#' @param interval Survival only: name of the aggregate interval-index
#'   column (values `1..K`) required when `cut_points` are supplied.
#' @param baseline Survival baseline hazard: `"piecewise"` (default,
#'   step function) or `"mspline"` (smooth M-spline over the `cut_points`
#'   grid).
#' @param n_basis Number of M-spline basis functions when
#'   `baseline = "mspline"` (must be `<=` the number of intervals).
#' @param cor Optional covariate correlation matrix for the Gaussian-copula
#'   integration. Defaults to the IPD correlation; pass a matrix to override
#'   or the identity to integrate covariates independently.
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
  if (as.integer(n_int) < 1L) {
    stop("`n_int` must be a positive integer.", call. = FALSE)
  }
  miss_em <- setdiff(effect_modifiers, names(ipd))
  if (length(miss_em)) {
    stop("`ipd` is missing effect-modifier column(s): ",
         paste(miss_em, collapse = ", "), call. = FALSE)
  }
  for (v in effect_modifiers) {
    if (!paste0(v, "_mean") %in% names(agd)) {
      stop("`agd` is missing `", v, "_mean`.", call. = FALSE)
    }
    if (family != "gaussian") {
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
  }

  # Reference treatment for reported relative effects.
  reference <- if (!is.null(inactive) && inactive %in% rownames(C)) {
    inactive
  } else {
    rownames(C)[1]
  }

  # Identifiability: component effects are estimable only if the within-study
  # arm contrasts span all components; otherwise some are prior-driven.
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
  rk <- if (is.null(Dmat) || !nrow(Dmat)) 0L else
    as.integer(Matrix::rankMatrix(Dmat))
  if (rk < ncol(C)) {
    warning("Only ", rk, " of ", ncol(C), " component effects are identified ",
            "by within-study contrasts; the remainder are informed only by ",
            "the prior. Provide evidence that contrasts those components.",
            call. = FALSE)
  }

  # Covariate correlation for the Gaussian-copula integration, estimated
  # from the (pre-expansion) IPD unless supplied via `cor`.
  cor_mat <- .cpaic_copula_cor(
    as.matrix(ipd[, effect_modifiers, drop = FALSE]), cor)

  # Survival: Lexis-expand IPD into interval-at-risk rows and configure the
  # piecewise-exponential baseline (K = 1 reduces to the exponential model).
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
      .cpaic_integration_points(means, sds, n_int_eff, cor = cor_mat)
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
    seed = seed %||% 1L, refresh = 0, show_messages = FALSE)
  .cpaic_check_diagnostics(fit)

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
