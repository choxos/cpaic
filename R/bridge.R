# Component-NMA bridge: reconnect a disconnected network ----------------------

#' Reconnect a network through its additive component structure
#'
#' Fits the additive component network meta-analysis (cNMA) model of Rücker
#' et al. (2020) to the aggregate contrast data, using
#' [netmeta::discomb()]. When the network is disconnected but its
#' sub-networks share components, the additive model estimates component
#' effects and so reconstructs relative effects *across* sub-networks. This
#' is the "connect first" step; population adjustment is layered on by
#' [cmaic()] / [cstc()], which replace unadjusted contrasts with adjusted
#' ones before calling this function.
#'
#' Estimability is checked per contrast, not by a single global rank test: a
#' relative effect is uniquely estimable if and only if its contrast vector
#' lies in the row space of the component design matrix `X = B C` (Wigle et
#' al. 2026). A rank-deficient network is therefore *not* rejected outright;
#' the contrasts that remain estimable are still reported, and those that are
#' not are returned as `NA` rather than as pseudoinverse artefacts. See
#' [estimable_effects()].
#'
#' @param network A [cpaic_network()] object.
#' @param common,random Fit common- and/or random-effects models.
#' @param ... Additional arguments passed to [netmeta::discomb()] (e.g.
#'   `tau.preset`).
#'
#' @return An object of class `cpaic_bridge` wrapping the
#'   [netmeta::discomb()] fit, with tidied component and treatment effects.
#' @references
#' Rücker G, Petropoulou M, Schwarzer G (2020). Network meta-analysis of
#' multicomponent interventions. *Biometrical Journal*, 62(3), 808--821.
#'
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
#' br <- cnma_bridge(net)
#' component_effects(br)
#' relative_effects(br)
#' @export
cnma_bridge <- function(network, common = FALSE, random = TRUE, ...) {
  stopifnot(inherits(network, "cpaic_network"))
  cols <- network$cols
  agd <- network$agd

  if (any(!is.finite(agd[[cols$TE]])) ||
      any(!is.finite(agd[[cols$seTE]]) | agd[[cols$seTE]] <= 0)) {
    stop("All contrasts must have a finite `TE` and a finite positive ",
         "`seTE` before bridging; NA/Inf/non-positive values remain ",
         "(an unadjusted IPD edge, perhaps). Fill them with cmaic()/cstc() ",
         "or supply complete aggregate data.", call. = FALSE)
  }

  conn <- cpaic_connectivity(network)
  has_components <- any(grepl(network$sep.comps, network$treatments,
                              fixed = TRUE))
  est <- conn$estimable
  if (!any(est$estimable)) {
    stop("No relative effect versus the reference (\"", network$reference,
         "\") is estimable from this network: every contrast lies outside ",
         "the row space of the component design matrix X = B C (rank(X) = ",
         conn$rank, " < ", conn$n_components, " components). The ",
         "sub-networks do not share enough components to bridge the gap. ",
         "See cpaic_connectivity() and estimable_effects().", call. = FALSE)
  }
  if (!all(est$estimable)) {
    warning("Some relative effects are not uniquely estimable and will be ",
            "returned as NA: ",
            paste(est$treatment[!est$estimable], collapse = ", "),
            ". They lie outside the row space of X (rank(X) = ", conn$rank,
            " < ", conn$n_components, " components). See ",
            "estimable_effects().", call. = FALSE)
  } else if (!conn$identifiable && has_components) {
    # Every requested contrast is estimable, but individual component effects
    # are not separately identified; those components are reported as NA.
    warning("All relative effects versus the reference are estimable, but ",
            "some individual component effects are not separately identified ",
            "(rank(X) = ", conn$rank, " < ", conn$n_components,
            " components) and are returned as NA: ",
            paste(names(conn$estimable_components)[
              !conn$estimable_components], collapse = ", "),
            ".", call. = FALSE)
  }

  run <- function() {
    netmeta::discomb(
      TE = agd[[cols$TE]],
      seTE = agd[[cols$seTE]],
      treat1 = agd[[cols$treat1]],
      treat2 = agd[[cols$treat2]],
      studlab = agd[[cols$studlab]],
      sm = network$sm,
      inactive = network$inactive,
      sep.comps = network$sep.comps,
      reference.group = network$reference,
      common = common,
      random = random,
      ...
    )
  }
  # For a singleton-treatment network discomb's component machinery is
  # trivially degenerate and emits noise; suppress it there (the treatment
  # effects are still correct). Real component networks keep their warnings.
  fit <- if (has_components) run() else suppressWarnings(run())

  effect <- if (random) "random" else "common"
  comp_tbl <- component_table(fit, effect = effect)
  # Blank out component effects the design cannot identify.
  bad_comp <- names(conn$estimable_components)[!conn$estimable_components]
  if (length(bad_comp)) {
    idx <- comp_tbl$component %in% bad_comp
    comp_tbl[idx, c("estimate", "se", "lower", "upper", "statistic",
                    "pval")] <- NA_real_
  }

  structure(
    list(
      fit = fit,
      components = comp_tbl,
      C.matrix = fit$C.matrix,
      sm = network$sm,
      reference = network$reference,
      effect = effect,
      Q = .cpaic_q_table(fit),
      connectivity = conn,
      estimable = est,
      network = network
    ),
    class = "cpaic_bridge"
  )
}

#' Collect the additive-model fit statistics from a discomb fit
#'
#' `Q.additive` is the *total* lack of fit of the additive component model.
#' The nested test of the additivity restrictions themselves is
#' `Q.diff = Q.additive - Q.standard`, which exists only when a standard
#' (non-additive) NMA is estimable, i.e. on a connected network. Neither
#' statistic can test whether component effects are constant *across*
#' disconnected sub-networks: there is no cross-gap evidence to test against.
#' @noRd
.cpaic_q_table <- function(fit) {
  g <- function(nm) {
    v <- fit[[nm]]
    if (is.null(v) || !length(v)) NA_real_ else as.numeric(v)[1]
  }
  list(
    Q = g("Q.additive"), df = g("df.Q.additive"), pval = g("pval.Q.additive"),
    Q.standard = g("Q.standard"), df.standard = g("df.Q.standard"),
    Q.diff = g("Q.diff"), df.diff = g("df.Q.diff"),
    pval.diff = g("pval.Q.diff")
  )
}

#' Tidy component-effect table from a discomb fit
#' @noRd
component_table <- function(fit, effect = c("random", "common")) {
  effect <- match.arg(effect)
  suffix <- if (effect == "random") "random" else "fixed"
  pick <- function(stub) {
    v <- fit[[paste0(stub, ".", suffix)]]
    as.vector(v)
  }
  out <- data.frame(
    component = fit$comps,
    estimate = pick("Comp"),
    se = pick("seComp"),
    lower = pick("lower.Comp"),
    upper = pick("upper.Comp"),
    statistic = pick("statistic.Comp"),
    pval = pick("pval.Comp"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "backtransf") <- isTRUE(fit$backtransf)
  attr(out, "sm") <- fit$sm
  out
}

#' @export
print.cpaic_bridge <- function(x, ...) {
  cat("cpaic component-NMA bridge (", x$effect, " effects, sm = ", x$sm,
      ")\n", sep = "")
  q <- x$Q
  cat("  Additive-model fit (Cochran Q): Q = ", round(q$Q, 2), ", df = ",
      q$df, ", p = ", format.pval(q$pval, digits = 3), "\n", sep = "")
  if (is.finite(q$Q.diff)) {
    cat("  Additivity restrictions (Q.diff): Q = ", round(q$Q.diff, 2),
        ", df = ", q$df.diff, ", p = ", format.pval(q$pval.diff, digits = 3),
        "\n", sep = "")
  }
  est <- x$estimable
  if (!all(est$estimable)) {
    cat("  Not estimable (NA): ",
        paste(est$treatment[!est$estimable], collapse = ", "), "\n", sep = "")
  }
  cat("\nComponent effects (", x$sm, " scale, link/log):\n", sep = "")
  comp <- x$components
  comp[, c("estimate", "se", "lower", "upper")] <-
    round(comp[, c("estimate", "se", "lower", "upper")], 3)
  print(comp[, c("component", "estimate", "se", "lower", "upper", "pval")],
        row.names = FALSE)
  invisible(x)
}

#' Component effects from a cpaic fit
#'
#' @param object A fitted cpaic object (`cpaic_bridge`, `cpaic_maic`,
#'   `cpaic_stc`, or `cpaic_mlnmr`).
#' @param newdata For [cmlnmr()] fits: a one-row data frame giving the target
#'   population's effect-modifier values, at which the component effects
#'   `beta + Gamma x` are reported. With `newdata = NULL` (default) the
#'   component *main* effects `beta` are returned; these are the effects at the
#'   covariate origin and are not population-adjusted.
#' @param ... Passed to methods.
#' @return A data frame of component effects (estimate, se, CI, p-value).
#'   Components that the design cannot identify are returned as `NA`.
#' @export
component_effects <- function(object, newdata = NULL, ...) {
  UseMethod("component_effects")
}

#' @export
component_effects.cpaic_bridge <- function(object, newdata = NULL, ...) {
  object$components
}
