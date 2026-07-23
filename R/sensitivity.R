# Bridge fragility: sensitivity of a cross-gap contrast to un-testable drift ---

#' Bridge fragility statistics from a draw vector of the contrast
#'
#' `D` is a vector of posterior draws of the link-scale contrast, `loading` is
#' the L1 loading `sum |m_c|`, `threshold` the decision boundary, and
#' `plausible_drift` an optional per-component drift bound for the robustness
#' probability. Returns the per-draw fragility threshold and summaries. This is
#' the pure numeric core, kept separate so it can be checked without a fit.
#' @noRd
.cpaic_bft_stats <- function(D, loading, threshold = 0,
                             plausible_drift = NULL) {
  bft <- abs(D - threshold) / loading            # per-draw drift to the boundary
  p_robust <- if (is.null(plausible_drift)) NA_real_ else
    mean(abs(D - threshold) > plausible_drift * loading)
  list(
    contrast_median = stats::median(D),
    loading = loading,
    bft_median = stats::median(bft),
    bft_lower = unname(stats::quantile(bft, 0.025)),
    bft_upper = unname(stats::quantile(bft, 0.975)),
    plausible_drift = plausible_drift %||% NA_real_,
    p_robust = p_robust,
    bft = bft
  )
}

#' Bridge fragility: how much cross-sub-network drift would change a conclusion
#'
#' In a disconnected network the cross-gap contrast exists only because the
#' component effects are assumed constant *across* sub-networks. That assumption
#' cannot be tested from the data, because there is no cross-gap evidence.
#' `bridge_fragility()` quantifies how sensitive a requested contrast is to a
#' violation of it.
#'
#' On the linear-predictor scale the contrast is \eqn{D = m'(\beta + \Gamma x)}
#' with \eqn{m = C_t - C_u}. A cross-sub-network drift \eqn{\Delta} in the
#' component effects shifts it to \eqn{D + m'\Delta}. Bounding each component's
#' drift by \eqn{|\Delta_c| \le d}, the worst-case shift is \eqn{d \sum_c |m_c|},
#' so the smallest per-component drift that moves the contrast to a decision
#' threshold \eqn{\tau} (default 0, on the link scale) is the **bridge fragility
#' threshold**
#' \deqn{\mathrm{BFT} = |D - \tau| / \textstyle\sum_c |m_c|,}
#' reported per posterior draw. A small BFT means a clinically trivial amount of
#' un-testable drift would overturn the conclusion. This is a conservative
#' worst-case over the component *main-effect* drift; interaction drift
#' \eqn{\Lambda} is not included, so the true fragility is no larger than
#' reported.
#'
#' @param object A [cmlnmr()] fit.
#' @param treatment,comparator The contrast to assess. `comparator` defaults to
#'   the fit reference.
#' @param newdata A one-row data frame giving the target population's
#'   effect-modifier values (required when the model has effect modifiers).
#' @param threshold Decision boundary on the link scale. Default `0` (no effect).
#' @param plausible_drift Optional per-component drift bound (link scale) at
#'   which to report the posterior probability that the conclusion is robust.
#' @param ... Unused.
#'
#' @return An object of class `cpaic_fragility`: the contrast, the L1 drift
#'   loading, the posterior of the bridge fragility threshold, and (if
#'   `plausible_drift` is given) the probability the conclusion survives it.
#' @seealso [relative_effects()], [estimable_effects_at()]
#' @examplesIf FALSE
#' bridge_fragility(fit, treatment = "A+B", newdata = data.frame(x1 = 0))
#' @export
bridge_fragility <- function(object, treatment, comparator = NULL,
                             newdata = NULL, threshold = 0,
                             plausible_drift = NULL, ...) {
  stopifnot(inherits(object, "cpaic_mlnmr"))
  C <- object$C.matrix
  trts <- rownames(C)
  if (is.null(comparator)) comparator <- object$reference
  if (!treatment %in% trts || !comparator %in% trts) {
    stop("`treatment` and `comparator` must be network treatments: ",
         paste(trts, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold)) {
    stop("`threshold` must be a finite number (on the link scale).",
         call. = FALSE)
  }
  if (!is.null(plausible_drift) &&
      (!is.numeric(plausible_drift) || length(plausible_drift) != 1L ||
       !is.finite(plausible_drift) || plausible_drift <= 0)) {
    stop("`plausible_drift` must be a positive number or NULL.", call. = FALSE)
  }
  x <- .cpaic_target_x(newdata, object$effect_modifiers)

  m <- C[treatment, ] - C[comparator, ]
  # Estimability at this target: a non-estimable contrast has no meaningful
  # fragility (the number would be the prior).
  v <- .cpaic_target_vec(m, x)
  N <- .cpaic_null_space(object$joint_design)
  if (!.cpaic_in_rowspace(matrix(v, nrow = 1L), N)) {
    stop("The contrast ", treatment, " vs ", comparator, " is not estimable at ",
         "this target population, so its bridge fragility is undefined. See ",
         "estimable_effects_at().", call. = FALSE)
  }

  loading <- sum(abs(m))
  if (loading == 0) {
    stop("The contrast has zero component loading (identical treatments); ",
         "bridge fragility is undefined.", call. = FALSE)
  }
  Beff <- .cpaic_beta_at(object, x)         # draws x components of (beta + Gamma x)
  D <- as.numeric(Beff %*% m)               # link-scale contrast draws

  stats <- .cpaic_bft_stats(D, loading, threshold, plausible_drift)
  out <- structure(
    c(list(treatment = treatment, comparator = comparator,
           threshold = threshold, sm = object$sm,
           target = stats::setNames(x, object$effect_modifiers)),
      stats[setdiff(names(stats), "bft")],
      list(bft_draws = stats$bft, contrast_draws = D)),
    class = "cpaic_fragility")
  out
}

#' @export
print.cpaic_fragility <- function(x, digits = 3, ...) {
  cat("Bridge fragility (", x$treatment, " vs ", x$comparator, ")\n", sep = "")
  if (length(x$target)) {
    cat("  Target population: ",
        paste(names(x$target), signif(x$target, 3), sep = " = ",
              collapse = ", "), "\n", sep = "")
  }
  cat("  Link-scale contrast (median): ", round(x$contrast_median, digits),
      " (threshold ", x$threshold, ")\n", sep = "")
  cat("  Worst-case drift loading (sum |m_c|): ", round(x$loading, digits),
      "\n", sep = "")
  cat("  Bridge fragility threshold (per-component drift to the boundary):\n")
  cat("    median ", round(x$bft_median, digits), " [", round(x$bft_lower, digits),
      ", ", round(x$bft_upper, digits), "] link-scale units\n", sep = "")
  if (is.finite(x$p_robust)) {
    cat("  P(conclusion robust to drift <= ", x$plausible_drift, "): ",
        round(x$p_robust, digits), "\n", sep = "")
  }
  cat("  Small values mean a clinically trivial, un-testable cross-sub-network\n",
      "  drift would overturn the conclusion.\n", sep = "")
  invisible(x)
}
