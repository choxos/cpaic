# Plots for the Bayesian component-additive ML-NMR ----------------------------
#
# Rankograms, deviance and leverage plots, prior-versus-posterior overlays,
# integration-error plots, MCMC diagnostics, and survival curves are ported from
# multinma (Phillippo et al. 2020, <doi:10.1111/rssa.12579>): specifically
# plot.nma_summary(), plot.nma_rank_probs(), plot.nma_dic(),
# plot_prior_posterior(), plot_integration_error(), plot.stan_nma(),
# plot.surv_nma_summary(), and geom_km(). The logic is re-implemented here on
# ggplot2 alone (multinma builds on ggdist and ggraph, which cpaic does not
# depend on). Both packages are licensed under GPL-3.
#
# The population-dependent rank curve and the estimability map have no
# counterpart in multinma: under population adjustment both the hierarchy and
# the estimable set are functions of the target population.

utils::globalVariables(c("dens", "density", "group", "parameter", "point",
                         "study_arm", "xseq"))

#' Posterior draws of one or more model parameters, as a plain matrix
#'
#' The result is a base matrix, not a `posterior::draws_matrix`. That matters:
#' `[.draws_matrix` defaults to `drop = FALSE`, so `m[i, j]` on a draws matrix
#' silently stays two-dimensional and downstream arithmetic goes wrong.
#'
#' Errors are not swallowed: a missing variable or an unreadable draws file is a
#' real problem, and reporting it as "nothing to plot" would hide it.
#' @noRd
.cpaic_draws <- function(object, variables) {
  out <- tryCatch(
    object$fit$draws(variables, format = "draws_matrix"),
    error = function(e) {
      stop("Could not read posterior draws for '",
           paste(variables, collapse = ", "), "' from the fit: ",
           conditionMessage(e),
           "\nA cmdstanr fit keeps its draws in temporary CSV files. If the fit ",
           "was saved and reloaded, save it with fit$fit$save_object() (or call ",
           "fit$fit$draws() before saving) so the draws travel with it.",
           call. = FALSE)
    })
  out <- as.matrix(out)
  res <- matrix(as.numeric(out), nrow = nrow(out), ncol = ncol(out))
  colnames(res) <- colnames(out)
  res
}

#' Check that we were handed a cmlnmr() fit
#' @noRd
.cpaic_check_mlnmr <- function(object, what) {
  if (!inherits(object, "cpaic_mlnmr")) {
    stop(what, " needs a cmlnmr() fit.", call. = FALSE)
  }
  invisible(TRUE)
}

# Rank probabilities: rankogram and cumulative ranks ---------------------------

#' Posterior rank probabilities in a target population
#'
#' The full rank distribution behind [cpaic_ranks()]: the posterior probability
#' that each treatment (or component) takes each rank, **in a named target
#' population**. Ported from `multinma::posterior_rank_probs()` (Phillippo et
#' al. 2020) and extended, because under population adjustment the hierarchy is
#' a function of the target: the component effects are `beta + Gamma x`, so the
#' ranks move with `x`.
#'
#' Elements that are not estimable at the target population are dropped from the
#' ranking set rather than ranked from the prior, exactly as in [cpaic_ranks()]
#' (Step 3 of Wigle et al. 2026); they are listed in the `dropped` attribute.
#'
#' @param object A [cmlnmr()] fit.
#' @param newdata A one-row data frame giving the target population's
#'   effect-modifier values.
#' @param what `"treatment"` (default) or `"component"`.
#' @param lower_is_better If `TRUE`, a smaller effect is preferred.
#' @param cumulative Return cumulative rank probabilities (the quantity SUCRA
#'   summarizes) instead of the rankogram? Default `FALSE`.
#' @param ... Unused.
#'
#' @return A data frame of class `cpaic_rank_probs` with one row per (element,
#'   rank) and columns `element`, `rank_position`, and `probability`.
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @seealso [cpaic_ranks()], [rank_curve()], [plot.cpaic_rank_probs()]
#' @examplesIf FALSE
#' rp <- rank_probs(fit, newdata = data.frame(x1 = 0.5), what = "component")
#' plot(rp)
#' @export
rank_probs <- function(object, newdata = NULL,
                       what = c("treatment", "component"),
                       lower_is_better = FALSE, cumulative = FALSE, ...) {
  .cpaic_check_mlnmr(object, "rank_probs()")
  what <- match.arg(what)
  if (!is.logical(cumulative) || length(cumulative) != 1L ||
      is.na(cumulative)) {
    stop("`cumulative` must be TRUE or FALSE.", call. = FALSE)
  }
  C <- object$C.matrix
  x <- .cpaic_target_x(newdata, object$effect_modifiers)
  Beff <- .cpaic_beta_at(object, x)

  if (what == "treatment") {
    ref <- object$reference
    elems <- setdiff(rownames(C), ref)
    Draws <- Beff %*% t(C)
    colnames(Draws) <- rownames(C)
    Draws <- Draws[, elems, drop = FALSE] - Draws[, ref]
    Lmat <- C[elems, , drop = FALSE] -
      matrix(C[ref, ], nrow = length(elems), ncol = ncol(C), byrow = TRUE)
  } else {
    elems <- object$comps
    Draws <- Beff
    colnames(Draws) <- elems
    Lmat <- diag(ncol(C))
    rownames(Lmat) <- elems
  }

  # Only rank what the target population actually identifies.
  V <- do.call(rbind, lapply(seq_along(elems),
                             function(i) .cpaic_target_vec(Lmat[i, ], x)))
  ok <- .cpaic_in_rowspace(V, .cpaic_null_space(object$joint_design))
  dropped <- elems[!ok]
  if (length(dropped)) {
    warning("Dropped from the hierarchy as not estimable in this target ",
            "population: ", paste(dropped, collapse = ", "),
            ". Ranking them would rank the prior. See estimable_effects_at().",
            call. = FALSE)
  }
  elems <- elems[ok]
  if (length(elems) < 2L) {
    stop("Fewer than two elements are estimable in this target population, so ",
         "no hierarchy can be formed. See estimable_effects_at().",
         call. = FALSE)
  }
  Draws <- Draws[, elems, drop = FALSE]

  sgn <- if (lower_is_better) 1 else -1
  R <- t(apply(sgn * Draws, 1L, rank, ties.method = "min"))
  n <- length(elems)
  P <- vapply(seq_len(n), function(k) colMeans(R == k), numeric(n))
  dimnames(P) <- list(elems, as.character(seq_len(n)))
  if (cumulative) P <- t(apply(P, 1L, cumsum))

  out <- data.frame(
    element = rep(elems, times = n),
    rank_position = rep(seq_len(n), each = n),
    probability = as.numeric(P),
    row.names = NULL, stringsAsFactors = FALSE)
  attr(out, "dropped") <- dropped
  attr(out, "target") <- x
  attr(out, "what") <- what
  attr(out, "cumulative") <- cumulative
  class(out) <- c("cpaic_rank_probs", "data.frame")
  out
}

#' Rankogram and cumulative rank plot
#'
#' Ported from `multinma::plot.nma_rank_probs()` (Phillippo et al. 2020). The
#' rankogram gives the posterior probability of each rank; the cumulative
#' version gives the probability of being ranked among the best `k`, whose
#' normalized area is SUCRA.
#'
#' Both are computed **in a named target population**, because a
#' population-adjusted hierarchy is not population-free.
#'
#' @param x A `cpaic_rank_probs` object from [rank_probs()].
#' @param y Unused, for compatibility with the [plot()] generic.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @seealso [rank_probs()], [plot_rank_curve()]
#' @examplesIf FALSE
#' plot(rank_probs(fit, newdata = data.frame(x1 = 0), what = "component"))
#' @export
plot.cpaic_rank_probs <- function(x, y, ...) {
  .cpaic_need_ggplot("plot.cpaic_rank_probs()")
  cumulative <- isTRUE(attr(x, "cumulative"))
  target <- attr(x, "target")
  dat <- as.data.frame(x)
  dat$element <- factor(dat$element, levels = unique(dat$element))
  nrank <- max(dat$rank_position)

  subtitle <- if (length(target)) {
    paste0("Target population: ",
           paste(names(target), signif(target, 3), sep = " = ",
                 collapse = ", "))
  } else {
    NULL
  }
  dropped <- attr(x, "dropped")
  caption <- if (length(dropped)) {
    paste0("Not estimable in this population, so not ranked: ",
           paste(dropped, collapse = ", "), ".")
  } else {
    NULL
  }

  p <- ggplot2::ggplot(
    dat, ggplot2::aes(x = rank_position, y = probability))

  p <- if (cumulative) {
    p + ggplot2::geom_line(colour = "#113259", linewidth = 0.7) +
      ggplot2::geom_point(colour = "#113259", size = 1.6)
  } else {
    p + ggplot2::geom_col(fill = "#55A480", colour = "grey30",
                          linewidth = 0.2, width = 0.8)
  }

  p +
    ggplot2::facet_wrap(~ element) +
    ggplot2::scale_x_continuous(
      "Rank", breaks = seq_len(nrank), minor_breaks = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      y = if (cumulative) "Cumulative probability" else "Probability",
      subtitle = subtitle, caption = caption) +
    .cpaic_theme()
}

#' Plot a population-adjusted hierarchy
#'
#' Plots the ranking metrics of a [cpaic_ranks()] hierarchy. Ranking metrics
#' depend on the set being ranked, so read them alongside the relative effects,
#' never instead of them.
#'
#' @param x A `cpaic_ranks` object from [cpaic_ranks()].
#' @param y Unused, for compatibility with the [plot()] generic.
#' @param ... Unused.
#' @param metric Which metric to plot: `"sucra"` (default), `"mean_rank"`,
#'   `"median_rank"`, or `"p_best"`.
#' @return A `ggplot` object.
#' @seealso [cpaic_ranks()], [plot_rank_curve()]
#' @examplesIf FALSE
#' plot(cpaic_ranks(fit, newdata = data.frame(x1 = 0)))
#' @export
plot.cpaic_ranks <- function(x, y, ...,
                             metric = c("sucra", "mean_rank", "median_rank",
                                        "p_best")) {
  .cpaic_need_ggplot("plot.cpaic_ranks()")
  metric <- match.arg(metric)
  target <- attr(x, "target")
  dropped <- attr(x, "dropped")
  dat <- as.data.frame(x)
  dat$value <- dat[[metric]]
  # Better is to the right: high SUCRA / p_best, low mean or median rank.
  decreasing <- metric %in% c("mean_rank", "median_rank")
  dat <- dat[order(dat$value, decreasing = decreasing), , drop = FALSE]
  dat$element <- factor(dat$element, levels = unique(dat$element))

  ggplot2::ggplot(dat, ggplot2::aes(x = value, y = element)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = value, y = element, yend = element),
      colour = "grey70") +
    ggplot2::geom_point(size = 3, colour = "#113259") +
    ggplot2::labs(
      x = switch(metric, sucra = "SUCRA", mean_rank = "Mean rank",
                 median_rank = "Median rank", p_best = "P(best)"),
      y = attr(x, "what"),
      subtitle = if (length(target)) {
        paste0("Target population: ",
               paste(names(target), signif(target, 3), sep = " = ",
                     collapse = ", "))
      },
      caption = if (length(dropped)) {
        paste0("Not estimable in this population, so not ranked: ",
               paste(dropped, collapse = ", "), ".")
      }) +
    .cpaic_theme()
}

# The headline cpaic plot: a hierarchy that moves with the target population ---

#' How the hierarchy changes across target populations
#'
#' **The headline figure of cpaic.** Under population adjustment the component
#' effects are `beta + Gamma x`, so a component's rank is a function of the
#' target population `x` and components **cross**: the component that leads in
#' one population can trail in another. A single hierarchy, quoted without a
#' population, is therefore not a well-posed answer. This plot shows the whole
#' family of hierarchies at once.
#'
#' There is no counterpart in multinma, which ranks in one population at a time.
#'
#' @param x A [cmlnmr()] fit, or the data frame returned by [rank_curve()].
#' @param em Name of the effect modifier to vary. Required when `x` is a fit.
#' @param values Numeric vector of target values for `em`. Required when `x` is
#'   a fit.
#' @param at Optional named vector fixing the other effect modifiers.
#' @param what,lower_is_better See [cpaic_ranks()].
#' @param metric Which ranking metric to trace: `"sucra"` (default),
#'   `"mean_rank"`, or `"p_best"`.
#' @param ... Passed to [rank_curve()] when `x` is a fit.
#'
#' @return A `ggplot` object.
#' @seealso [rank_curve()], [cpaic_ranks()], [plot_estimability()]
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @examplesIf FALSE
#' plot_rank_curve(fit, em = "x1", values = seq(-1, 1, by = 0.25),
#'                 what = "component")
#' @export
plot_rank_curve <- function(x, em = NULL, values = NULL, at = NULL,
                            what = c("treatment", "component"),
                            lower_is_better = FALSE,
                            metric = c("sucra", "mean_rank", "p_best"), ...) {
  .cpaic_need_ggplot("plot_rank_curve()")
  metric <- match.arg(metric)

  if (inherits(x, "cpaic_mlnmr")) {
    what <- match.arg(what)
    if (is.null(em) || is.null(values)) {
      stop("`em` and `values` are required when `x` is a cmlnmr() fit: name ",
           "the effect modifier to vary and the target values to vary it over.",
           call. = FALSE)
    }
    curve <- rank_curve(x, em = em, values = values, at = at, what = what,
                        lower_is_better = lower_is_better, ...)
  } else if (is.data.frame(x) && !is.null(attr(x, "em"))) {
    curve <- x
  } else {
    stop("`x` must be a cmlnmr() fit or the data frame returned by ",
         "rank_curve().", call. = FALSE)
  }

  em <- attr(curve, "em")
  what <- attr(curve, "what")
  dat <- as.data.frame(curve)
  dat$target <- dat[[em]]
  dat$value <- dat[[metric]]
  dat$element <- factor(dat$element, levels = unique(dat$element))

  ylab <- switch(metric, sucra = "SUCRA", mean_rank = "Mean rank",
                 p_best = "P(best)")

  ggplot2::ggplot(dat, ggplot2::aes(x = target, y = value, colour = element)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_viridis_d(
      paste0(toupper(substring(what, 1, 1)), substring(what, 2)),
      end = 0.9) +
    ggplot2::labs(
      x = paste0("Target population (", em, ")"),
      y = ylab,
      subtitle = paste0("The hierarchy is a function of the target ",
                        "population; where the lines cross, the ",
                        what, " ordering reverses"),
      caption = paste0("Ranking metrics depend on the set ranked. Report ",
                       "them with the relative effects, not instead of them.")
    ) +
    .cpaic_theme()
}

# Estimability across target populations ---------------------------------------

#' Map which contrasts are estimable, and on what evidence, across populations
#'
#' Estimability under population adjustment is itself a function of the target
#' population: a contrast identified at the covariate origin need not be
#' identified in a population where the relevant component by effect-modifier
#' interactions are not pinned down. This plot evaluates
#' [estimable_effects_at()] over a grid of target populations and tiles the
#' result, separating contrasts identified by **IPD** (a within-study
#' interaction, which randomization protects) from those identified only
#' **ecologically**, from between-study differences in aggregate covariate means
#' (which randomization does not protect; Berlin et al. 2002).
#'
#' There is no counterpart in multinma.
#'
#' @param object A [cmlnmr()] fit.
#' @param em Name of the effect modifier to vary across the grid.
#' @param values Numeric vector of target values for `em`.
#' @param at Optional named vector fixing the other effect modifiers. Defaults
#'   to 0 for each.
#' @param reference Reference treatment. Defaults to the fit's reference.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @seealso [estimable_effects_at()], [plot_rank_curve()]
#' @references
#' Berlin JA, Santanna J, Schmid CH, Szczech LA, Feldman HI (2002).
#' Individual patient- versus group-level data meta-regressions for the
#' investigation of treatment effect modifiers. *Statistics in Medicine*,
#' 21(3), 371--387.
#' @examplesIf FALSE
#' plot_estimability(fit, em = "x1", values = seq(-1, 1, by = 0.5))
#' @export
plot_estimability <- function(object, em, values, at = NULL, reference = NULL,
                              ...) {
  .cpaic_need_ggplot("plot_estimability()")
  .cpaic_check_mlnmr(object, "plot_estimability()")
  ems <- object$effect_modifiers
  if (!is.character(em) || length(em) != 1L || !em %in% ems) {
    stop("`em` must be one of the effect modifiers: ",
         paste(ems, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(values) || !length(values) || any(!is.finite(values))) {
    stop("`values` must be a non-empty vector of finite target values.",
         call. = FALSE)
  }
  base <- stats::setNames(rep(0, length(ems)), ems)
  if (!is.null(at)) base[names(at)] <- at

  res <- lapply(values, function(v) {
    nd <- as.data.frame(as.list(replace(base, em, v)))
    tab <- estimable_effects_at(object, newdata = nd, reference = reference)
    tab$target <- v
    tab
  })
  dat <- do.call(rbind, res)
  dat$contrast <- paste0(dat$treatment, " vs ", dat$comparator)
  dat$contrast <- factor(dat$contrast, levels = rev(unique(dat$contrast)))
  dat$identified_by <- factor(
    dat$identified_by, levels = c("IPD", "aggregate", "none"),
    labels = c("IPD (within-study)", "Aggregate (ecological)",
               "Not estimable"))

  ggplot2::ggplot(
    dat, ggplot2::aes(x = target, y = contrast, fill = identified_by)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(
      "Identified by",
      values = c(`IPD (within-study)` = "#55A480",
                 `Aggregate (ecological)` = "#E8C468",
                 `Not estimable` = "#B2182B"),
      drop = FALSE) +
    ggplot2::labs(
      x = paste0("Target population (", em, ")"),
      y = "Contrast",
      subtitle = "Estimability is a function of the target population",
      caption = paste0("An ecologically identified contrast rests on ",
                       "between-study covariate means, which randomization ",
                       "does not protect.")) +
    .cpaic_theme()
}

# Deviance, dev-dev, and leverage ---------------------------------------------

#' Residual deviance draws, fitted values, and observations, per data point
#'
#' The saturated log likelihood is subtracted so that a well-fitting data point
#' contributes about one residual deviance, as in the DIC diagnostics of
#' multinma and NICE TSD2. For the binomial and Poisson families it is a
#' per-point constant, so it cancels from DIC *differences* and leaves the `pV`
#' penalty of [dic()] untouched. For the Gaussian family the IPD saturated term
#' depends on `sigma`, so it is a draw-level quantity.
#' @noRd
.cpaic_resdev_draws <- function(object) {
  .cpaic_check_mlnmr(object, "This diagnostic")
  args <- object$refit_args
  if (is.null(args)) {
    stop("The fit does not carry the data needed for deviance diagnostics.",
         call. = FALSE)
  }
  ll <- .cpaic_draws(object, "log_lik")
  if (!ncol(ll)) {
    stop("The fit does not contain pointwise `log_lik` draws.", call. = FALSE)
  }
  ipd <- args$ipd
  agd <- args$agd
  n_ipd <- nrow(ipd)
  n_agd <- nrow(agd)
  fam <- object$family

  labels <- c(
    if (n_ipd) paste0(as.character(ipd[[args$study]]), ": ",
                      as.character(ipd[[args$trt]]), " [", seq_len(n_ipd), "]"),
    if (n_agd) paste0(as.character(agd[[args$study]]), ": ",
                      as.character(agd[[args$trt]]))
  )
  type <- c(rep("IPD", n_ipd), rep("AgD", n_agd))

  # Saturated log likelihood: the best any model could do for this data point.
  # For a Bernoulli IPD point the saturated fit is exact, so its saturated log
  # likelihood is zero.
  D <- -2 * ll                                   # deviance draws
  saturated <- TRUE
  if (fam == "binomial") {
    r <- as.numeric(agd[[args$r]]); n <- as.numeric(agd[[args$n]])
    sat <- c(rep(0, n_ipd), stats::dbinom(r, n, r / n, log = TRUE))
    D <- sweep(D, 2L, 2 * sat, "+")
  } else if (fam == "poisson") {
    yy <- as.numeric(ipd[[args$outcome]])
    r <- as.numeric(agd[[args$r]])
    sat <- c(stats::dpois(yy, pmax(yy, .Machine$double.eps), log = TRUE),
             stats::dpois(r, pmax(r, .Machine$double.eps), log = TRUE))
    D <- sweep(D, 2L, 2 * sat, "+")
  } else if (fam == "gaussian") {
    # The IPD saturated log likelihood depends on sigma, so it is a draw-level
    # quantity rather than a per-point constant.
    sig <- as.numeric(.cpaic_draws(object, "sigma"))
    if (n_ipd) {
      D[, seq_len(n_ipd)] <- D[, seq_len(n_ipd), drop = FALSE] +
        2 * (-log(sig * sqrt(2 * pi)))
    }
    if (n_agd) {
      se <- as.numeric(agd[[args$se]])
      j <- n_ipd + seq_len(n_agd)
      D[, j] <- sweep(D[, j, drop = FALSE], 2L,
                      2 * stats::dnorm(0, 0, se, log = TRUE), "+")
    }
  } else {
    saturated <- FALSE                           # survival: no saturated model
  }

  list(deviance = D, labels = labels, type = type, saturated = saturated,
       family = fam, n_ipd = n_ipd, n_agd = n_agd)
}

#' Posterior mean fitted value and observed value, per data point
#' @noRd
.cpaic_fitted_observed <- function(object) {
  args <- object$refit_args
  ipd <- args$ipd
  agd <- args$agd
  n_ipd <- nrow(ipd); n_agd <- nrow(agd)
  fam <- object$family
  mn <- function(v) colMeans(.cpaic_draws(object, v))
  switch(
    fam,
    binomial = list(
      observed = c(as.numeric(ipd[[args$outcome]]), as.numeric(agd[[args$r]])),
      fitted = c(stats::plogis(mn("eta_ipd")),
                 mn("p_agd") * as.numeric(agd[[args$n]]))),
    gaussian = list(
      observed = c(as.numeric(ipd[[args$outcome]]),
                   as.numeric(agd[[args$outcome]])),
      fitted = c(mn("eta_ipd"), mn("eta_agd"))),
    poisson = list(
      observed = c(as.numeric(ipd[[args$outcome]]), as.numeric(agd[[args$r]])),
      fitted = c(exp(mn("eta_ipd")), mn("lambda_agd"))),
    list(observed = rep(NA_real_, n_ipd + n_agd),
         fitted = rep(NA_real_, n_ipd + n_agd))
  )
}

#' Deviance and dev-dev plots
#'
#' Ported from `multinma::plot.nma_dic()` (Phillippo et al. 2020). With a single
#' [dic()] object the plot shows each data point's contribution to the posterior
#' mean deviance; points contributing much more than the rest are fitted poorly.
#' With two [dic()] objects it draws the **dev-dev plot**: points below the line
#' of equality fit better under the second model, points above it fit better
#' under the first.
#'
#' [dic()] stores the posterior mean deviance per data point, not the residual
#' deviance. The two differ by the saturated log likelihood, which for the
#' binomial and Poisson families is a function of the data alone: it is the same
#' under both models, so it shifts both axes equally and the line of equality
#' keeps its meaning. For a Gaussian likelihood the saturated term also involves
#' the model's own `sigma`, so two Gaussian models with very different residual
#' variance are shifted by different amounts; read that comparison with care.
#'
#' For posterior uncertainty and for the leverage plot, which need the saturated
#' model explicitly, call [plot_leverage()] on the fitted model itself.
#'
#' @param x A `cpaic_dic` object from [dic()].
#' @param y An optional second `cpaic_dic` object, for a dev-dev plot.
#' @param ... Unused.
#' @param labels Model names used on the axes of a dev-dev plot.
#'
#' @return A `ggplot` object.
#' @seealso [dic()], [plot_leverage()]
#' @examplesIf FALSE
#' plot(dic(fit_fe), dic(fit_re), labels = c("Fixed", "Random"))
#' @export
plot.cpaic_dic <- function(x, y = NULL, ..., labels = c("Model 1", "Model 2")) {
  .cpaic_need_ggplot("plot.cpaic_dic()")
  if (!inherits(x, "cpaic_dic")) {
    stop("`x` must be a cpaic_dic object from dic().", call. = FALSE)
  }

  if (is.null(y)) {
    dat <- data.frame(point = seq_along(x$pointwise),
                      deviance_x = as.numeric(x$pointwise))
    return(
      ggplot2::ggplot(dat, ggplot2::aes(x = point, y = deviance_x)) +
        ggplot2::geom_point(size = 1.8, colour = "#113259") +
        ggplot2::labs(
          x = "Data point", y = "Posterior mean deviance contribution",
          subtitle = paste0("DIC = ", round(x$dic, 1), " (pV = ",
                            round(x$p_eff, 1), ")")) +
        .cpaic_theme()
    )
  }

  if (!inherits(y, "cpaic_dic")) {
    stop("`y` must be a cpaic_dic object from dic(), or NULL.", call. = FALSE)
  }
  if (length(x$pointwise) != length(y$pointwise)) {
    stop("The two models were fitted to different numbers of data points, so ",
         "their deviance contributions cannot be compared point by point.",
         call. = FALSE)
  }
  if (!is.character(labels) || length(labels) != 2L) {
    stop("`labels` must be a character vector of length 2.", call. = FALSE)
  }

  dat <- data.frame(deviance_x = as.numeric(x$pointwise),
                    deviance_y = as.numeric(y$pointwise))
  ulim <- max(dat$deviance_x, dat$deviance_y, na.rm = TRUE)
  llim <- min(0, dat$deviance_x, dat$deviance_y, na.rm = TRUE)

  ggplot2::ggplot(dat, ggplot2::aes(x = deviance_x, y = deviance_y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    ggplot2::geom_point(size = 1.8, alpha = 0.75, colour = "#113259") +
    ggplot2::coord_fixed(xlim = c(llim, ulim), ylim = c(llim, ulim)) +
    ggplot2::labs(
      x = paste0("Deviance contribution (", labels[1], ")"),
      y = paste0("Deviance contribution (", labels[2], ")"),
      caption = paste0("Below the line: better fit under ", labels[2],
                       ". Above: better fit under ", labels[1], ".")) +
    .cpaic_theme()
}

#' Leverage plot
#'
#' Ported from `multinma::plot.nma_dic()` with `type = "leverage"` (Phillippo et
#' al. 2020). Each data point's leverage (its contribution to the effective
#' number of parameters) is plotted against its signed square root residual
#' deviance, with contours of constant DIC contribution. Points outside the
#' `DIC = 3` contour are usually the ones spoiling the fit.
#'
#' Leverage is the pointwise `pV` penalty, half the posterior variance of the
#' point's deviance, matching the penalty used by [dic()]. Because deviances
#' covary across points, the pointwise leverages need not sum exactly to the
#' model's total `pV`.
#'
#' Leverage plots need a saturated model and so are **not available for survival
#' outcomes**, where the censored contributions have no saturated reference.
#' multinma declines them for the same reason.
#'
#' @param object A [cmlnmr()] fit.
#' @param ... Unused.
#' @param dic_contours Numeric vector of DIC contours to draw. Default `1:4`.
#'
#' @return A `ggplot` object.
#' @seealso [dic()], [plot.cpaic_dic()]
#' @examplesIf FALSE
#' plot_leverage(fit)
#' @export
plot_leverage <- function(object, ..., dic_contours = 1:4) {
  .cpaic_need_ggplot("plot_leverage()")
  rd <- .cpaic_resdev_draws(object)
  if (!rd$saturated) {
    stop("Leverage plots need a saturated model, which censored survival data ",
         "do not have. Use plot.cpaic_dic() to compare deviance contributions ",
         "between models instead.", call. = FALSE)
  }
  if (!is.numeric(dic_contours) || any(!is.finite(dic_contours))) {
    stop("`dic_contours` must be a numeric vector.", call. = FALSE)
  }

  fo <- .cpaic_fitted_observed(object)
  resdev <- colMeans(rd$deviance)
  leverage <- apply(rd$deviance, 2L, stats::var) / 2
  sgn <- sign(fo$observed - fo$fitted)
  sgn[sgn == 0] <- 1

  dat <- data.frame(
    .label = rd$labels, type = rd$type,
    resdev = resdev, leverage = leverage,
    ssrd = sgn * sqrt(pmax(resdev, 0)),
    stringsAsFactors = FALSE)
  rmax <- max(abs(dat$ssrd), sqrt(max(dic_contours)), na.rm = TRUE) * 1.05

  contours <- do.call(rbind, lapply(dic_contours, function(cc) {
    xs <- seq(-rmax, rmax, length.out = 201)
    data.frame(ssrd = xs, leverage = cc - xs^2, group = paste0("DIC = ", cc),
               stringsAsFactors = FALSE)
  }))

  ggplot2::ggplot(dat, ggplot2::aes(x = ssrd, y = leverage)) +
    ggplot2::geom_line(
      data = contours,
      mapping = ggplot2::aes(x = ssrd, y = leverage, group = group),
      colour = "grey60", linewidth = 0.3, inherit.aes = FALSE) +
    ggplot2::geom_point(ggplot2::aes(colour = type), size = 2) +
    ggplot2::scale_colour_manual(
      "Evidence", values = c(AgD = "#113259", IPD = "#55A480"), drop = FALSE) +
    ggplot2::coord_cartesian(
      xlim = c(-rmax, rmax),
      ylim = c(min(0, dat$leverage, na.rm = TRUE),
               max(dic_contours, dat$leverage, na.rm = TRUE))) +
    ggplot2::labs(x = "Signed square root residual deviance", y = "Leverage",
                  caption = paste0("Contours of constant DIC contribution: ",
                                   paste(dic_contours, collapse = ", "), ".")) +
    .cpaic_theme()
}

# Prior versus posterior --------------------------------------------------------

#' Density of a cpaic prior on a grid
#' @noRd
.cpaic_prior_density <- function(spec, n = 501L) {
  dist <- spec$distribution
  scale <- spec$scale
  df <- spec$df
  positive <- dist %in% c("half-normal", "half-student-t")
  q <- switch(
    dist,
    normal = stats::qnorm(c(0.001, 0.999), 0, scale),
    student_t = scale * stats::qt(c(0.001, 0.999), df),
    `half-normal` = c(0, stats::qnorm(0.999, 0, scale)),
    `half-student-t` = c(0, scale * stats::qt(0.999, df)),
    stop("Unsupported prior distribution '", dist, "'.", call. = FALSE)
  )
  xseq <- seq(q[1], q[2], length.out = n)
  dens <- switch(
    dist,
    normal = stats::dnorm(xseq, 0, scale),
    student_t = stats::dt(xseq / scale, df) / scale,
    `half-normal` = 2 * stats::dnorm(xseq, 0, scale),
    `half-student-t` = 2 * stats::dt(xseq / scale, df) / scale
  )
  if (positive) dens[xseq < 0] <- 0
  data.frame(xseq = xseq, dens = dens, stringsAsFactors = FALSE)
}

#' Prior versus posterior
#'
#' Ported from `multinma::plot_prior_posterior()` (Phillippo et al. 2020).
#' Posteriors are drawn as histograms, priors as lines. Where a posterior simply
#' reproduces its prior, the data carry no information about that parameter, and
#' any quantity that leans on it is prior-driven rather than estimated. This is
#' the visual counterpart of [prior_sensitivity()].
#'
#' It matters most for the component by effect-modifier interactions `gamma`:
#' interactions informed only by aggregate arms are weakly identified, and
#' `prior_gamma_scale` then does real regularization.
#'
#' @param x A [cmlnmr()] fit.
#' @param ... Unused.
#' @param prior Which priors to show. Any of `"intercept"` (`mu`), `"beta"`
#'   (component effects), `"regression"` (`breg`), `"gamma"` (component by
#'   effect-modifier interactions), and `"tau"` (heterogeneity, random effects
#'   only). Defaults to all that the model used.
#' @param bins Number of histogram bins for the posterior. Default `40`.
#'
#' @return A `ggplot` object.
#' @seealso [prior_sensitivity()], [prior_predictive_check()]
#' @examplesIf FALSE
#' plot_prior_posterior(fit, prior = "gamma")
#' @export
plot_prior_posterior <- function(x, ..., prior = NULL, bins = 40) {
  .cpaic_need_ggplot("plot_prior_posterior()")
  .cpaic_check_mlnmr(x, "plot_prior_posterior()")
  par_of <- c(intercept = "mu", beta = "beta", regression = "breg",
              gamma = "gamma", tau = "tau")
  used <- names(par_of)
  if (identical(x$trt_effects, "fixed")) used <- setdiff(used, "tau")
  if (is.null(prior)) {
    prior <- used
  } else if (!is.character(prior) || !all(prior %in% used)) {
    stop("`prior` must be a character vector with elements from: ",
         paste(used, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(bins) || length(bins) != 1L || bins < 2) {
    stop("`bins` must be a single number of at least 2.", call. = FALSE)
  }

  post <- list()
  pri <- list()
  for (nm in prior) {
    spec <- x$priors[[nm]]
    d <- .cpaic_draws(x, par_of[[nm]])
    if (!ncol(d)) next
    dens <- .cpaic_prior_density(spec)
    for (j in seq_len(ncol(d))) {
      pname <- colnames(d)[j]
      post[[length(post) + 1L]] <- data.frame(
        parameter = pname, value = as.numeric(d[, j]),
        stringsAsFactors = FALSE)
      pri[[length(pri) + 1L]] <- data.frame(
        parameter = pname, xseq = dens$xseq, dens = dens$dens,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(post)) {
    stop("No prior/posterior pair to plot.", call. = FALSE)
  }
  post <- do.call(rbind, post)
  pri <- do.call(rbind, pri)
  lev <- unique(post$parameter)
  post$parameter <- factor(post$parameter, levels = lev)
  pri$parameter <- factor(pri$parameter, levels = lev)

  ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data = post,
      ggplot2::aes(x = value, y = ggplot2::after_stat(density)),
      bins = bins, fill = "grey75", colour = "grey45", linewidth = 0.15) +
    ggplot2::geom_line(
      data = pri, ggplot2::aes(x = xseq, y = dens),
      colour = "#B2182B", linewidth = 0.7) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::labs(
      x = "Value", y = "Density",
      subtitle = "Posterior (histogram) against prior (red line)",
      caption = paste0("A posterior that reproduces its prior is not ",
                       "estimated from the data. See prior_sensitivity().")) +
    .cpaic_theme()
}

# Integration error -------------------------------------------------------------

#' Numerical integration error against the number of integration points
#'
#' Ported from `multinma::plot_integration_error()` (Phillippo et al. 2020).
#' Aggregate arms are fitted by integrating the individual-level model over the
#' study's covariate distribution with Sobol' quasi-Monte-Carlo points. The
#' integration error at `N` points is the estimate using the first `N` points
#' minus the estimate using all of them; the typical convergence rate of QMC
#' integration, `1/N`, is drawn for reference. If the error has not settled well
#' inside the `1/N` envelope by `n_int`, refit with more integration points.
#'
#' cpaic does not save cumulative integration points inside Stan, so the
#' integrated aggregate-arm quantity is reconstructed here from the posterior
#' draws and the (deterministic) Sobol' point set. This is exact, not an
#' approximation, but it is not free: subsample the posterior with `ndraws` on a
#' large fit.
#'
#' Not available for `family = "survival"`, where the aggregate contribution is
#' a `log_sum_exp` over integration points of the likelihood rather than an
#' integrated mean outcome; multinma declines this plot for survival models too.
#' Nor for a Gaussian model fitted with all-normal margins, which is *exact* at
#' the covariate means and uses a single integration point.
#'
#' @param x A [cmlnmr()] fit.
#' @param ... Unused.
#' @param int_thin Report the error every `int_thin` points. Default is
#'   `n_int / 8`, rounded up.
#' @param ndraws Number of posterior draws to summarize over. Default `200`.
#' @param show_expected_rate Draw the `1/N` convergence envelope? Default
#'   `TRUE`.
#'
#' @return A `ggplot` object.
#' @seealso [cmlnmr()]
#' @examplesIf FALSE
#' plot_integration_error(fit)
#' @export
plot_integration_error <- function(x, ..., int_thin = NULL, ndraws = 200L,
                                   show_expected_rate = TRUE) {
  .cpaic_need_ggplot("plot_integration_error()")
  .cpaic_check_mlnmr(x, "plot_integration_error()")
  if (identical(x$family, "survival")) {
    stop("Integration-error plots are not supported for survival models: the ",
         "aggregate contribution is a log_sum_exp of the likelihood over the ",
         "integration points, not an integrated mean outcome.", call. = FALSE)
  }
  args <- x$refit_args
  agd <- args$agd
  ems <- x$effect_modifiers
  margins <- x$margins
  n_int <- as.integer(args$n_int)
  if (identical(x$family, "gaussian") && all(margins == "normal")) {
    stop("A Gaussian model with normal margins is integrated exactly at the ",
         "covariate means (one integration point), so there is no integration ",
         "error to plot.", call. = FALSE)
  }
  if (n_int < 2L) {
    stop("The model was fitted with a single integration point; there is no ",
         "integration error to trace.", call. = FALSE)
  }

  # Rebuild the (deterministic) Sobol' integration points used at fit time.
  sd_of <- function(v, a) {
    if (margins[[v]] == "bernoulli") return(NA_real_)
    agd[[paste0(v, "_sd")]][a]
  }
  X_list <- lapply(seq_len(nrow(agd)), function(a) {
    means <- vapply(ems, function(v) agd[[paste0(v, "_mean")]][a], numeric(1))
    sds <- vapply(ems, sd_of, numeric(1), a = a)
    .cpaic_integration_points(means, sds, n_int, cor = x$cor,
                              margins = margins)
  })

  mu <- .cpaic_draws(x, "mu")
  beta <- .cpaic_draws(x, "beta")
  breg <- .cpaic_draws(x, "breg")
  gamma <- .cpaic_draws(x, "gamma")
  nd <- min(as.integer(ndraws), nrow(beta))
  keep <- unique(round(seq(1, nrow(beta), length.out = nd)))
  studies <- sort(unique(c(as.character(args$ipd[[args$study]]),
                           as.character(agd[[args$study]]))))
  Tc <- x$C.matrix[match(as.character(agd[[args$trt]]),
                         rownames(x$C.matrix)), , drop = FALSE]
  nC <- ncol(x$C.matrix)
  Q <- length(ems)
  invlink <- switch(x$family, binomial = stats::plogis, poisson = exp,
                    gaussian = identity)

  if (is.null(int_thin)) int_thin <- max(1L, ceiling(n_int / 8))
  int_thin <- as.integer(int_thin)
  if (int_thin < 1L || int_thin > n_int) {
    stop("`int_thin` must be between 1 and the number of integration points.",
         call. = FALSE)
  }
  grid <- unique(c(seq(int_thin, n_int, by = int_thin), n_int))

  re <- .cpaic_random_effects(args$ipd, agd, args$study, args$trt)
  delta <- if (identical(x$trt_effects, "random")) {
    .cpaic_draws(x, "delta")
  } else {
    NULL
  }

  out <- list()
  for (a in seq_len(nrow(agd))) {
    Xa <- X_list[[a]]
    s <- match(as.character(agd[[args$study]][a]), studies)
    # Linear predictor at every integration point, for every kept draw.
    tc <- Tc[a, ]
    # gamma is stored as gamma[component, effect modifier].
    gcols <- matrix(0, nrow = length(keep), ncol = Q)
    for (q in seq_len(Q)) {
      cols <- paste0("gamma[", seq_len(nC), ",", q, "]")
      gcols[, q] <- as.numeric(gamma[keep, cols, drop = FALSE] %*% tc)
    }
    base <- mu[keep, s] + as.numeric(beta[keep, , drop = FALSE] %*% tc)
    if (!is.null(delta) && re$re_idx_agd[a] > 0L) {
      base <- base + delta[keep, re$re_idx_agd[a]]
    }
    # eta[draw, point] = base + X %*% breg + X %*% gcols (row-wise).
    eta <- outer(base, rep(1, n_int)) +
      breg[keep, , drop = FALSE] %*% t(Xa) +
      gcols %*% t(Xa)
    mean_out <- invlink(eta)                     # natural scale
    cum <- t(apply(mean_out, 1L, cumsum)) /
      matrix(seq_len(n_int), nrow = length(keep), ncol = n_int, byrow = TRUE)
    final <- cum[, n_int]
    for (g in grid[grid < n_int]) {
      out[[length(out) + 1L]] <- data.frame(
        study_arm = paste0(agd[[args$study]][a], ": ", agd[[args$trt]][a]),
        n_int = g, value = cum[, g] - final,
        stringsAsFactors = FALSE)
    }
  }
  dat <- do.call(rbind, out)
  if (is.null(dat) || !nrow(dat)) {
    stop("No integration error to plot; increase `n_int` or lower `int_thin`.",
         call. = FALSE)
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = n_int, y = value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60")

  if (isTRUE(show_expected_rate)) {
    conv <- data.frame(
      n_int = rep(seq(1, n_int, length.out = 201), 2L),
      value = c(seq(1, n_int, length.out = 201)^-1,
                -seq(1, n_int, length.out = 201)^-1),
      group = rep(c("pos", "neg"), each = 201L),
      stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_line(
      data = conv, ggplot2::aes(x = n_int, y = value, group = group),
      colour = "grey60", linetype = 2, inherit.aes = FALSE)
  }

  p +
    ggplot2::geom_point(alpha = 0.15, size = 0.7, colour = "#113259") +
    ggplot2::coord_cartesian(ylim = range(dat$value, na.rm = TRUE)) +
    ggplot2::facet_wrap(~ study_arm) +
    ggplot2::labs(x = "Number of integration points",
                  y = "Estimated integration error",
                  caption = "Dashed lines: the typical 1/N QMC rate.") +
    .cpaic_theme()
}

# MCMC diagnostics --------------------------------------------------------------

#' MCMC diagnostics for a cML-NMR fit
#'
#' Ported from `multinma::plot.stan_nma()` (Phillippo et al. 2020), which hands
#' the posterior draws to `bayesplot`. Diverging transitions, a high maximum
#' `Rhat`, or a low effective sample size all mean the posterior has not been
#' explored, and nothing downstream of it is trustworthy.
#'
#' @param x A [cmlnmr()] fit.
#' @param y Unused, for compatibility with the [plot()] generic.
#' @param ... Passed to the underlying `bayesplot` function.
#' @param type `"trace"` (default), `"density"`, `"hist"`, `"pairs"`, `"rhat"`,
#'   or `"neff"`.
#' @param pars Character vector of parameter names to show. Defaults to the
#'   component effects `beta`, the interactions `gamma`, and, for a
#'   random-effects model, the heterogeneity `tau`.
#'
#' @return A `ggplot` object (`"pairs"` returns a `bayesplot` grid).
#' @seealso [cmlnmr()], [plot_prior_posterior()]
#' @examplesIf FALSE
#' plot(fit, type = "trace")
#' plot(fit, type = "pairs", pars = c("beta[1]", "tau[1]"))
#' @export
plot.cpaic_mlnmr <- function(x, y, ...,
                             type = c("trace", "density", "hist", "pairs",
                                      "rhat", "neff"),
                             pars = NULL) {
  .cpaic_need_ggplot("plot.cpaic_mlnmr()")
  type <- match.arg(type)
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("MCMC diagnostic plots need the 'bayesplot' package. Install it ",
         "with install.packages(\"bayesplot\").", call. = FALSE)
  }
  if (is.null(pars)) {
    pars <- c("beta", "gamma")
    if (identical(x$trt_effects, "random")) pars <- c(pars, "tau")
  }
  if (!is.character(pars) || !length(pars)) {
    stop("`pars` must be a non-empty character vector.", call. = FALSE)
  }

  if (type %in% c("rhat", "neff")) {
    s <- x$fit$summary(pars)
    if (type == "rhat") {
      return(bayesplot::mcmc_rhat(stats::setNames(s$rhat, s$variable), ...) +
               .cpaic_theme())
    }
    ndraws <- .cpaic_ndraws(x)
    return(bayesplot::mcmc_neff(
      stats::setNames(s$ess_bulk / ndraws, s$variable), ...) + .cpaic_theme())
  }

  draws <- x$fit$draws(pars, format = "draws_array")
  fn <- switch(type,
               trace = bayesplot::mcmc_trace,
               density = bayesplot::mcmc_dens_overlay,
               hist = bayesplot::mcmc_hist,
               pairs = bayesplot::mcmc_pairs)
  p <- fn(draws, ...)
  if (type == "pairs") return(p)
  p + .cpaic_theme()
}

#' Total number of posterior draws in a cmdstanr fit
#' @noRd
.cpaic_ndraws <- function(object) nrow(.cpaic_draws(object, "beta"))

# Survival curves ----------------------------------------------------------------

#' Kaplan-Meier curves from the survival data behind a cML-NMR fit
#'
#' Ported from `multinma::geom_km()` (Phillippo et al. 2020). Returns a list of
#' ggplot2 layers, so it can be added to an existing plot; adding it to
#' [plot_survival()] overlays the observed data on the fitted survival curves.
#'
#' Only status `1` counts as an event; statuses `0`, `2`, and `3` (right, left,
#' and interval censoring) are treated as censored for the empirical curve.
#'
#' @param object A [cmlnmr()] fit with `family = "survival"`.
#' @param ... Passed to [survival::survfit()].
#' @param curve_args,cens_args Optional lists of arguments customizing the
#'   curves ([ggplot2::geom_step()]) and the censoring marks
#'   ([ggplot2::geom_point()]).
#'
#' @return A list of ggplot2 layers.
#' @seealso [plot_survival()]
#' @examplesIf FALSE
#' plot_survival(fit) + geom_km(fit)
#' @export
geom_km <- function(object, ..., curve_args = list(), cens_args = list()) {
  .cpaic_need_ggplot("geom_km()")
  .cpaic_check_mlnmr(object, "geom_km()")
  if (!identical(object$family, "survival")) {
    stop("geom_km() needs a cmlnmr() fit with family = \"survival\".",
         call. = FALSE)
  }
  dat <- .cpaic_surv_rows(object)

  km <- do.call(rbind, lapply(split(dat, dat$study_arm), function(d) {
    sf <- survival::survfit(
      survival::Surv(d$time, as.integer(d$status == 1L)) ~ 1, ...)
    rbind(
      data.frame(study_arm = d$study_arm[1], time = 0, surv = 1, n.censor = 0,
                 stringsAsFactors = FALSE),
      data.frame(study_arm = d$study_arm[1], time = sf$time, surv = sf$surv,
                 n.censor = sf$n.censor, stringsAsFactors = FALSE))
  }))
  km$study_arm <- factor(km$study_arm, levels = sort(unique(km$study_arm)))

  curve <- c(list(mapping = ggplot2::aes(x = time, y = surv,
                                         group = study_arm),
                  data = km, linewidth = 0.3, colour = "grey25"),
             curve_args)
  cens <- c(list(mapping = ggplot2::aes(x = time, y = surv,
                                        group = study_arm),
                 data = km[km$n.censor >= 1L, , drop = FALSE],
                 shape = 3, size = 1.2, stroke = 0.3, colour = "grey25"),
            cens_args)
  curve <- curve[!duplicated(names(curve))]
  cens <- cens[!duplicated(names(cens))]

  list(do.call(ggplot2::geom_step, curve),
       do.call(ggplot2::geom_point, cens))
}

#' Survival rows (IPD plus reconstructed aggregate rows) behind a fit
#' @noRd
.cpaic_surv_rows <- function(object) {
  args <- object$refit_args
  grab <- function(d) {
    if (is.null(d) || !nrow(d)) return(NULL)
    data.frame(
      study = as.character(d[[args$study]]),
      trt = as.character(d[[args$trt]]),
      time = as.numeric(d[[args$time]]),
      status = as.integer(d[[args$outcome]]),
      stringsAsFactors = FALSE)
  }
  out <- rbind(grab(args$ipd), grab(args$agd))
  out$study_arm <- paste0(out$study, ": ", out$trt)
  out
}

#' Fitted survival curves from a cML-NMR fit
#'
#' Ported from `multinma::plot.surv_nma_summary()` (Phillippo et al. 2020).
#' Draws the model-implied survival function for each study arm, averaged over
#' that arm's own covariate distribution (the IPD covariates for an IPD arm, the
#' integration points for an aggregate arm), with posterior credible bands. Add
#' [geom_km()] to overlay the observed Kaplan-Meier curves; systematic
#' departures indicate that the baseline hazard is too rigid.
#'
#' The baseline hazard is whatever the model was fitted with: a piecewise
#' exponential step function, or a continuous cubic M-spline. Its posterior
#' enters the curve through the integrated basis, so the bands include baseline
#' uncertainty.
#'
#' @param object A [cmlnmr()] fit with `family = "survival"`.
#' @param ... Unused.
#' @param times Time grid. Defaults to 100 points spanning the observed times.
#' @param ndraws Number of posterior draws to summarize over. Default `200`.
#' @param level Credible level for the band. Default `0.95`.
#'
#' @return A `ggplot` object.
#' @seealso [geom_km()], [cmlnmr()]
#' @examplesIf FALSE
#' plot_survival(fit) + geom_km(fit)
#' @export
plot_survival <- function(object, ..., times = NULL, ndraws = 200L,
                          level = 0.95) {
  .cpaic_need_ggplot("plot_survival()")
  .cpaic_check_mlnmr(object, "plot_survival()")
  if (!identical(object$family, "survival")) {
    stop("plot_survival() needs a cmlnmr() fit with family = \"survival\".",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  args <- object$refit_args
  rows <- .cpaic_surv_rows(object)
  if (is.null(times)) {
    times <- seq(min(rows$time) / 10, max(rows$time), length.out = 100L)
  }
  if (!is.numeric(times) || any(!is.finite(times)) || any(times <= 0)) {
    stop("`times` must be positive and finite.", call. = FALSE)
  }

  # Rebuild the baseline basis exactly as cmlnmr() did.
  basis_cuts <- if (identical(args$baseline, "mspline")) numeric() else
    args$cut_points
  spec <- .cpaic_survival_basis_spec(
    observed_times = rows$time, baseline = args$baseline,
    cut_points = basis_cuts, n_basis = args$n_basis)
  H0 <- .cpaic_survival_basis_eval(spec, times, rep(0, length(times)),
                                   rep(0, length(times)))$itime

  mu <- .cpaic_draws(object, "mu")
  beta <- .cpaic_draws(object, "beta")
  breg <- .cpaic_draws(object, "breg")
  gamma <- .cpaic_draws(object, "gamma")
  coefs <- .cpaic_draws(object, "coefficients")
  nd <- min(as.integer(ndraws), nrow(beta))
  keep <- unique(round(seq(1, nrow(beta), length.out = nd)))

  ems <- object$effect_modifiers
  Q <- length(ems)
  nC <- ncol(object$C.matrix)
  studies <- sort(unique(c(as.character(args$ipd[[args$study]]),
                           as.character(args$agd[[args$study]]))))
  delta <- if (identical(object$trt_effects, "random")) {
    .cpaic_draws(object, "delta")
  } else {
    NULL
  }

  # Covariate values supporting each arm.
  cov_of_arm <- .cpaic_arm_covariates(object, ems)
  arms <- names(cov_of_arm)
  a <- (1 - level) / 2

  res <- lapply(arms, function(key) {
    info <- cov_of_arm[[key]]
    Xa <- info$X
    tc <- object$C.matrix[info$trt, ]
    s <- match(info$study, studies)
    gcols <- matrix(0, nrow = length(keep), ncol = Q)
    for (q in seq_len(Q)) {
      cols <- paste0("gamma[", seq_len(nC), ",", q, "]")
      gcols[, q] <- as.numeric(gamma[keep, cols, drop = FALSE] %*% tc)
    }
    base <- mu[keep, s] + as.numeric(beta[keep, , drop = FALSE] %*% tc)
    if (!is.null(delta) && info$re_idx > 0L) {
      base <- base + delta[keep, info$re_idx]
    }
    eta <- outer(base, rep(1, nrow(Xa))) +
      breg[keep, , drop = FALSE] %*% t(Xa) + gcols %*% t(Xa)  # draws x n
    Hbase <- coefs[keep, , drop = FALSE] %*% t(H0)             # draws x times
    # S_bar(t) = mean_i exp(-H0(t) exp(eta_i)), averaged over the arm's
    # covariate distribution: the marginal survival the KM curve estimates.
    S <- vapply(seq_along(times), function(k) {
      rowMeans(exp(-Hbase[, k] * exp(eta)))
    }, numeric(length(keep)))
    q <- apply(S, 2L, stats::quantile, probs = c(a, 0.5, 1 - a), names = FALSE)
    data.frame(study_arm = key, time = times, surv = q[2, ],
               lower = q[1, ], upper = q[3, ], stringsAsFactors = FALSE)
  })
  dat <- do.call(rbind, res)
  dat$study_arm <- factor(dat$study_arm, levels = sort(unique(dat$study_arm)))

  ggplot2::ggplot(dat, ggplot2::aes(x = time, y = surv)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, group = study_arm),
      fill = "#55A480", alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(group = study_arm),
                       colour = "#113259", linewidth = 0.7) +
    ggplot2::facet_wrap(~ study_arm) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Time", y = "Survival probability",
      subtitle = paste0("Fitted survival (", args$baseline,
                        " baseline), averaged over each arm's covariates")) +
    .cpaic_theme()
}

#' Covariate values supporting each study arm, and its random-effect index
#' @noRd
.cpaic_arm_covariates <- function(object, ems) {
  args <- object$refit_args
  ipd <- args$ipd
  agd <- args$agd
  re <- .cpaic_random_effects(ipd, agd, args$study, args$trt)
  out <- list()

  if (!is.null(ipd) && nrow(ipd)) {
    key <- paste0(as.character(ipd[[args$study]]), ": ",
                  as.character(ipd[[args$trt]]))
    for (k in unique(key)) {
      idx <- which(key == k)
      out[[k]] <- list(
        study = as.character(ipd[[args$study]])[idx[1]],
        trt = as.character(ipd[[args$trt]])[idx[1]],
        X = as.matrix(ipd[idx, ems, drop = FALSE]),
        re_idx = re$re_idx_ipd[idx[1]])
    }
  }
  if (!is.null(agd) && nrow(agd)) {
    n_int <- as.integer(args$n_int)
    margins <- object$margins
    key <- paste0(as.character(agd[[args$study]]), ": ",
                  as.character(agd[[args$trt]]))
    sd_of <- function(v, a) {
      if (margins[[v]] == "bernoulli") return(NA_real_)
      agd[[paste0(v, "_sd")]][a]
    }
    for (k in unique(key)) {
      a1 <- which(key == k)[1]
      means <- vapply(ems, function(v) agd[[paste0(v, "_mean")]][a1],
                      numeric(1))
      sds <- vapply(ems, sd_of, numeric(1), a = a1)
      X <- .cpaic_integration_points(means, sds, n_int, cor = object$cor,
                                     margins = margins)
      colnames(X) <- ems
      out[[k]] <- list(
        study = as.character(agd[[args$study]])[a1],
        trt = as.character(agd[[args$trt]])[a1],
        X = X,
        re_idx = re$re_idx_agd[a1])
    }
  }
  out
}
