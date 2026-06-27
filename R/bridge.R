# Component-NMA bridge: reconnect a disconnected network ----------------------

#' Reconnect a network through its additive component structure
#'
#' Fits the additive component network meta-analysis (cNMA) model of Rücker
#' et al. (2020) to the aggregate contrast data, using
#' [netmeta::discomb()]. When the network is disconnected but its
#' sub-networks share components, the additive model estimates all
#' component effects and so reconstructs the (otherwise unavailable)
#' relative effects *across* sub-networks. This is the "connect first"
#' step; population adjustment is layered on by [cmaic()] / [cstc()], which
#' replace unadjusted contrasts with adjusted ones before calling this
#' function.
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

  conn <- cpaic_connectivity(network)
  has_components <- any(grepl(network$sep.comps, network$treatments,
                             fixed = TRUE))
  if (!conn$identifiable) {
    if (!conn$connected) {
      # Disconnected AND non-identifiable: the cross-sub-network contrasts
      # are not estimable; refuse rather than return spurious finite effects.
      stop("The network is disconnected and cannot be bridged: the ",
           "sub-networks do not share enough components to identify the ",
           "component effects (rank(X) = ", conn$rank, " < ",
           conn$n_components, "). See cpaic_connectivity().", call. = FALSE)
    } else if (has_components) {
      # Connected but a component is not separately identifiable; the
      # treatment effects are still defined, the component will be NA. (A
      # connected singleton-treatment network is component-degenerate by
      # construction and needs no warning.)
      warning("Component effects are not uniquely identifiable: rank(X) = ",
              conn$rank, " < ", conn$n_components, " components. Some effects ",
              "will be NA. See cpaic_connectivity().", call. = FALSE)
    }
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

  structure(
    list(
      fit = fit,
      components = comp_tbl,
      C.matrix = fit$C.matrix,
      sm = network$sm,
      reference = network$reference,
      effect = effect,
      Q = list(Q = fit$Q.additive, df = fit$df.Q.additive,
               pval = fit$pval.Q.additive),
      connectivity = conn,
      network = network
    ),
    class = "cpaic_bridge"
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
  cat("  Additivity (Cochran Q): Q = ", round(x$Q$Q, 2), ", df = ", x$Q$df,
      ", p = ", format.pval(x$Q$pval, digits = 3), "\n", sep = "")
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
#' @param object A `cpaic_bridge`, `cpaic_maic`, or `cpaic_stc` object.
#' @param ... Unused.
#' @return A data frame of component effects (estimate, se, CI, p-value).
#' @export
component_effects <- function(object, ...) {
  UseMethod("component_effects")
}

#' @export
component_effects.cpaic_bridge <- function(object, ...) object$components
