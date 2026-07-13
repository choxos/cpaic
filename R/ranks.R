# Population-adjusted treatment and component hierarchies ----------------------
#
# Wigle et al. (2026) set out a workflow for building treatment and component
# hierarchies in an (aggregate-data) component NMA:
#
#   Step 1  state the idealized hierarchy question, i.e. the set S to be ranked
#           and the criterion (the ranking metric);
#   Step 2  determine which of the required relative effects are estimable;
#   Step 3  refine S to the estimable subset, or decline to rank;
#   Step 4  compute the ranking metrics and report them alongside the effects.
#
# Under population adjustment the component effects are population-specific,
#
#     theta_t(x) = C_t' (beta + Gamma x),
#
# so the RANKING IS TOO. A component can lead in one population and trail in
# another, and the hierarchy question changes from "which component is best?" to
# "which component is best IN THIS POPULATION?".
#
# Step 2 acquires a population-adjusted analogue as well: the estimable set is
# itself a function of the target population (see R/estimability.R), so the set
# that can legitimately be ranked can differ between target populations. This
# function performs Steps 2 to 4 automatically and reports what it dropped.

#' Population-adjusted treatment and component hierarchies
#'
#' Ranks treatments or components *in a named target population*, following the
#' workflow of Wigle et al. (2026) but with every quantity evaluated at the
#' target's effect-modifier values. Because the component effects are
#' population-specific under population adjustment, so is the hierarchy: a
#' component may rank first in one population and last in another.
#'
#' Elements whose relative effect is not estimable at that target population are
#' **dropped from the ranking set** rather than ranked from a prior-driven
#' posterior, and are reported in the `dropped` attribute. This is Step 3 of the
#' Wigle et al. workflow, and it matters more here than in the aggregate-data
#' case, because the estimable set depends on the target (see
#' [estimable_effects_at()]).
#'
#' Ranking metrics depend on the set being ranked, so they are not comparable
#' across different sets. Report them alongside the relative effects, never
#' instead of them.
#'
#' @param object A [cmlnmr()] fit.
#' @param newdata A one-row data frame giving the target population's
#'   effect-modifier values. Required when the model has effect modifiers.
#' @param what `"treatment"` (default) or `"component"`. Ranking components by
#'   their incremental effect is only meaningful in an additive model.
#' @param set Optional character vector restricting the elements to rank (the
#'   set `S` of Wigle et al.). Defaults to all treatments (or all components).
#' @param lower_is_better If `TRUE`, a smaller effect is preferred (e.g.
#'   mortality). Default `FALSE` (a larger effect is preferred).
#' @param ... Unused.
#'
#' @return A data frame, ordered from most to least preferred, with columns
#'   `element`, `estimate` (posterior mean of the relative effect versus the
#'   reference, on the link scale), `p_best`, `median_rank`, `mean_rank` and
#'   `sucra`. The `dropped` attribute lists elements excluded as not estimable
#'   in this target population.
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @seealso [estimable_effects_at()], [relative_effects()]
#' @examplesIf FALSE
#' # Which component is best for a patient population with x1 = 0.5?
#' cpaic_ranks(fit, newdata = data.frame(x1 = 0.5), what = "component")
#' @export
cpaic_ranks <- function(object, newdata = NULL,
                        what = c("treatment", "component"), set = NULL,
                        lower_is_better = FALSE, ...) {
  stopifnot(inherits(object, "cpaic_mlnmr"))
  what <- match.arg(what)
  C <- object$C.matrix
  ems <- object$effect_modifiers
  x <- .cpaic_target_x(newdata, ems)

  # Draws of the population-specific component effects beta + Gamma x.
  Beff <- .cpaic_beta_at(object, x)

  N <- .cpaic_null_space(object$joint_design)
  if (what == "treatment") {
    ref <- object$reference
    elems <- setdiff(rownames(C), ref)
    Draws <- Beff %*% t(C)
    colnames(Draws) <- rownames(C)
    Draws <- Draws[, elems, drop = FALSE] - Draws[, ref]   # vs the reference
    Lmat <- C[elems, , drop = FALSE] -
      matrix(C[ref, ], nrow = length(elems), ncol = ncol(C), byrow = TRUE)
  } else {
    elems <- object$comps
    Draws <- Beff                       # incremental effect of each component
    colnames(Draws) <- elems
    Lmat <- diag(ncol(C))
    rownames(Lmat) <- elems
  }

  # Step 2: which of these are estimable IN THIS TARGET POPULATION?
  V <- do.call(rbind, lapply(seq_along(elems), function(i) {
    .cpaic_target_vec(Lmat[i, ], x)
  }))
  ok <- .cpaic_in_rowspace(V, N)
  names(ok) <- elems

  if (!is.null(set)) {
    bad <- setdiff(set, elems)
    if (length(bad)) {
      stop("`set` names elements not in the network: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    keep0 <- elems %in% set
    elems <- elems[keep0]; ok <- ok[keep0]
    Draws <- Draws[, elems, drop = FALSE]
  }

  # Step 3: refine the set to the estimable elements, and say what was dropped.
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

  # Step 4: ranking metrics, computed within the refined set.
  sgn <- if (lower_is_better) 1 else -1     # rank() is ascending
  R <- t(apply(sgn * Draws, 1L, rank, ties.method = "average"))
  n <- length(elems)

  out <- data.frame(
    element     = elems,
    estimate    = colMeans(Draws),
    p_best      = colMeans(R == 1),
    median_rank = apply(R, 2L, stats::median),
    mean_rank   = colMeans(R),
    sucra       = (n - colMeans(R)) / (n - 1),
    row.names   = NULL, stringsAsFactors = FALSE
  )
  out <- out[order(out$mean_rank), ]
  rownames(out) <- NULL
  attr(out, "dropped") <- dropped
  attr(out, "target") <- x
  attr(out, "what") <- what
  class(out) <- c("cpaic_ranks", "data.frame")
  out
}

#' @export
print.cpaic_ranks <- function(x, digits = 3, ...) {
  tgt <- attr(x, "target")
  cat("Population-adjusted ", attr(x, "what"), " hierarchy\n", sep = "")
  if (length(tgt)) {
    cat("  Target population: ",
        paste(names(tgt), signif(tgt, 3), sep = " = ", collapse = ", "),
        "\n", sep = "")
  }
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE)
  dr <- attr(x, "dropped")
  if (length(dr)) {
    cat("  Not estimable in this population, so not ranked: ",
        paste(dr, collapse = ", "), "\n", sep = "")
  }
  cat("  Ranking metrics depend on the set ranked; report them with the",
      " effects, not instead.\n", sep = "")
  invisible(x)
}

#' How a hierarchy changes across target populations
#'
#' Recomputes [cpaic_ranks()] over a grid of target populations, so that the
#' population dependence of the hierarchy is visible. Under population adjustment
#' a component's rank is a function of the target, and this is the object that
#' shows it.
#'
#' @param object A [cmlnmr()] fit.
#' @param em Name of the effect modifier to vary.
#' @param values Numeric vector of target values for `em`.
#' @param at Optional named vector fixing the other effect modifiers. Defaults
#'   to 0 for each.
#' @param what,lower_is_better See [cpaic_ranks()].
#' @param ... Unused.
#' @return A data frame with one row per (element, target value), giving `sucra`,
#'   `mean_rank`, `p_best` and `estimate`, plus `estimable`.
#' @seealso [cpaic_ranks()]
#' @examplesIf FALSE
#' rank_curve(fit, em = "x1", values = seq(-1, 1, by = 0.25), what = "component")
#' @export
rank_curve <- function(object, em, values, at = NULL,
                       what = c("treatment", "component"),
                       lower_is_better = FALSE, ...) {
  stopifnot(inherits(object, "cpaic_mlnmr"))
  what <- match.arg(what)
  ems <- object$effect_modifiers
  if (!em %in% ems) {
    stop("`em` must be one of the effect modifiers: ",
         paste(ems, collapse = ", "), call. = FALSE)
  }
  base <- stats::setNames(rep(0, length(ems)), ems)
  if (!is.null(at)) base[names(at)] <- at

  res <- lapply(values, function(v) {
    nd <- as.data.frame(as.list(replace(base, em, v)))
    r <- tryCatch(
      suppressWarnings(cpaic_ranks(object, newdata = nd, what = what,
                                   lower_is_better = lower_is_better)),
      error = function(e) NULL)
    if (is.null(r)) return(NULL)
    data.frame(element = r$element, value = v, estimate = r$estimate,
               p_best = r$p_best, mean_rank = r$mean_rank, sucra = r$sucra,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) {
    stop("No target value in `values` yields two or more estimable elements.",
         call. = FALSE)
  }
  names(out)[names(out) == "value"] <- em
  attr(out, "em") <- em
  attr(out, "what") <- what
  out
}
