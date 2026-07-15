# Diagnostics: additivity, effective sample size ------------------------------

#' Fit statistics for the additive component model
#'
#' Returns the Cochran Q statistics from [netmeta::discomb()].
#'
#' Two different statistics are reported, and only one of them tests
#' additivity:
#'
#' * `Q` (`Q.additive`) is the **total lack of fit** of the additive component
#'   model, pooling ordinary heterogeneity/inconsistency with any failure of
#'   additivity. It is not a test of additivity.
#' * `Q.diff = Q.additive - Q.standard` is the **nested test of the additivity
#'   restrictions** themselves, and is the statistic to read. It exists only
#'   when a standard (non-additive) NMA is also estimable, i.e. on a connected
#'   network; on a disconnected network it is `NA`.
#'
#' **Neither statistic can test the assumption that actually bridges a
#' disconnected network**, namely that component effects (and, under
#' population adjustment, component x effect-modifier interactions) are
#' *constant across sub-networks*. There is by construction no cross-gap
#' evidence against which to test it. A large p-value here is therefore not a
#' licence to bridge; the assumption must be defended on clinical grounds
#' (Veroniki et al. 2026).
#'
#' @param object A `cpaic_bridge` or `cpaic_fit` object.
#' @param ... Unused.
#' @return A one-row data frame with `Q`, `df`, `pval` (total additive-model
#'   lack of fit) and `Q.diff`, `df.diff`, `pval.diff` (the nested additivity
#'   test, `NA` on a disconnected network).
#' @export
additivity_test <- function(object, ...) {
  UseMethod("additivity_test")
}

#' @export
additivity_test.cpaic_fit <- function(object, ...) {
  additivity_test(object$bridge, ...)
}

#' @export
additivity_test.cpaic_mlnmr <- function(object, ...) {
  stop("additivity_test() (Cochran Q) applies to the frequentist component ",
       "bridge, not to cmlnmr(). For the Bayesian model, compare an additive ",
       "fit against one with interaction terms using information criteria ",
       "such as LOO via loo::loo(fit) or WAIC via loo::waic(fit).",
       call. = FALSE)
}

#' @export
additivity_test.cpaic_bridge <- function(object, ...) {
  q <- object$Q
  out <- data.frame(Q = q$Q, df = q$df, pval = q$pval,
                    Q.diff = q$Q.diff, df.diff = q$df.diff,
                    pval.diff = q$pval.diff, row.names = NULL)
  attr(out, "connected") <- isTRUE(object$connectivity$connected)
  class(out) <- c("cpaic_additivity", "data.frame")
  out
}

#' @export
print.cpaic_additivity <- function(x, ...) {
  cat("Additive component model: fit statistics\n")
  if (isTRUE(x$df == 0)) {
    # A saturated model has no residual degrees of freedom, so Q is identically
    # zero by arithmetic. Printing "Q = 0" invites the reader to see a perfect
    # fit, when in truth the statistic has NO power and cannot detect anything.
    # In simulation, cSTC coverage fell to 0.50 under a large additivity
    # violation while this statistic reported Q = 0 with df = 0.
    cat("  !! SATURATED MODEL: 0 residual degrees of freedom.\n")
    cat("     Q = 0 here is arithmetic, NOT evidence of fit. This statistic has\n")
    cat("     no power and cannot detect a violation of additivity, however\n")
    cat("     large. Do not read it as reassurance.\n")
    invisible(x)
    return(invisible(x))
  }
  cat("  Total lack of fit (Q.additive): Q = ", round(x$Q, 3), ", df = ",
      x$df, ", p = ", format.pval(x$pval, digits = 3), "\n", sep = "")
  if (is.finite(x$Q.diff)) {
    cat("  Additivity restrictions (Q.diff = Q.additive - Q.standard):\n")
    cat("    Q = ", round(x$Q.diff, 3), ", df = ", x$df.diff, ", p = ",
        format.pval(x$pval.diff, digits = 3), "\n", sep = "")
    if (!is.na(x$pval.diff) && x$pval.diff < 0.05) {
      cat("    -> evidence against additivity within the network.\n")
    }
  } else {
    cat("  Additivity restrictions (Q.diff): not available; no standard NMA\n",
        "    is estimable on a disconnected network.\n", sep = "")
  }
  cat("  Note: neither statistic tests whether component effects are constant\n",
      "  ACROSS sub-networks, which is the assumption that bridges the gap.\n",
      "  That assumption is untestable from the data and must be justified\n",
      "  clinically.\n", sep = "")
  invisible(x)
}

#' Does the individual patient data actually inform this contrast?
#'
#' Population adjustment only helps if the adjusted edges actually *carry* the
#' contrast you care about. They need not. In a component bridge the estimate of
#' a contrast `m' beta` is a weighted combination of the observed edges,
#'
#' \deqn{m'\hat\beta = \underbrace{m' (X'WX)^{+} X'W}_{w} \, d ,}
#'
#' so edge `j` influences the answer only through its weight `w_j`. An IPD edge
#' with `w_j` of zero contributes nothing to that contrast, and adjusting it
#' changes nothing.
#'
#' The weight uses a diagonal `W` of inverse edge variances. The fit itself is
#' produced by `netmeta::discomb()`, which accounts for the within-study
#' covariance of a multi-arm trial, so in a network containing multi-arm studies
#' the weight reported here is a close approximation to the fitted estimator's
#' influence rather than its exact value. It is intended as a screening
#' diagnostic, to flag an IPD edge that carries little or no weight on the
#' contrast; read a weight near zero as "this edge barely matters here", not as
#' an exact sensitivity.
#'
#' This matters because the usual diagnostic cannot detect the problem. In
#' simulation, putting the IPD on an edge that does not bridge the gap left
#' cMAIC numerically identical to the unadjusted analysis (bias +0.374, coverage
#' 0.676) while [effective_sample_size()] happily reported an ESS of 999 out of
#' 1000. A healthy ESS says the *weights* are well behaved; it says nothing about
#' whether the reweighted edge is relevant to your estimand.
#'
#' @param object A `cpaic_bridge`, `cpaic_maic` or `cpaic_stc` object.
#' @param treatment,comparator The contrast of interest. `comparator` defaults
#'   to the network reference.
#' @param tol Influence weights below this (relative to the largest) are treated
#'   as zero.
#' @param ... Unused.
#'
#' @return A data frame with one row per edge: `studlab`, `treat1`, `treat2`,
#'   `has_ipd`, and `influence` (the weight `w_j`). Edges are ordered by
#'   absolute influence. A warning is issued if any IPD edge has no influence on
#'   the requested contrast.
#' @seealso [effective_sample_size()], [estimable_effects()]
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
#' br <- cnma_bridge(net)
#' edge_influence(br, treatment = "A+B+C")
#' @export
edge_influence <- function(object, treatment, comparator = NULL, tol = 1e-8,
                           ...) {
  if (inherits(object, "cpaic_fit")) object <- object$bridge
  if (!inherits(object, "cpaic_bridge")) {
    stop("edge_influence() needs a cpaic_bridge (or a cMAIC / cSTC fit).",
         call. = FALSE)
  }
  conn <- object$connectivity
  net <- object$network
  cols <- net$cols
  C <- conn$C
  trts <- rownames(C)
  if (is.null(comparator)) comparator <- object$reference
  if (!treatment %in% trts || !comparator %in% trts) {
    stop("`treatment` and `comparator` must be network treatments.",
         call. = FALSE)
  }

  X <- conn$X                                  # edges x components
  seTE <- net$agd[[cols$seTE]]
  tau2 <- if (identical(object$effect, "random")) {
    t2 <- object$fit$tau2
    if (is.null(t2) || !is.finite(t2)) 0 else t2
  } else 0
  w <- 1 / (seTE^2 + tau2)                     # inverse-variance weights
  XtWX <- t(X) %*% (w * X)
  m <- C[treatment, ] - C[comparator, ]
  infl <- as.numeric(m %*% MASS_ginv(XtWX) %*% t(X) %*% diag(w, length(w)))

  ipd_studies <- if (is.null(net$ipd_info)) character(0) else
    net$ipd_info$studies
  out <- data.frame(
    studlab = as.character(net$agd[[cols$studlab]]),
    treat1 = as.character(net$agd[[cols$treat1]]),
    treat2 = as.character(net$agd[[cols$treat2]]),
    has_ipd = as.character(net$agd[[cols$studlab]]) %in% ipd_studies,
    influence = infl,
    row.names = NULL, stringsAsFactors = FALSE
  )
  # `tol` is relative to the largest influence, as documented. Flooring the
  # scale at 1 (as an earlier version did) makes the cutoff absolute whenever
  # every influence is below 1, which fires a false "no influence" warning: three
  # equal edges of 1/3 have a documented cutoff of tol/3, not tol. When no edge
  # has any influence the relative comparison is undefined, so nothing is flagged.
  scale <- max(abs(out$influence))
  dead <- out$has_ipd & scale > 0 & abs(out$influence) < tol * scale
  if (any(dead)) {
    warning("Individual patient data on ",
            paste(unique(out$studlab[dead]), collapse = ", "),
            " have NO influence on ", treatment, " versus ", comparator,
            ". Adjusting those edges cannot change this contrast, whatever the ",
            "effective sample size says. See ?edge_influence.", call. = FALSE)
  }
  out <- out[order(-abs(out$influence)), ]
  rownames(out) <- NULL
  out
}

#' Moore-Penrose pseudo-inverse (avoids a MASS dependency)
#' @noRd
MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  keep <- s$d > max(tol * s$d[1], 0)
  if (!any(keep)) return(matrix(0, ncol(A), nrow(A)))
  s$v[, keep, drop = FALSE] %*%
    ((1 / s$d[keep]) * t(s$u[, keep, drop = FALSE]))
}

#' Effective sample sizes from a cMAIC fit
#'
#' The effective sample size summarizes the precision lost to reweighting. It is
#' **not** a validity diagnostic: a healthy ESS says the weights are well
#' behaved, not that the reweighted edge is relevant to your estimand. Use
#' [edge_influence()] to ask whether the IPD informs the contrast at all.
#'
#' @param object A `cpaic_maic` object.
#' @param ... Unused.
#' @return A named numeric vector of effective sample sizes per IPD study.
#' @seealso [edge_influence()]
#' @export
effective_sample_size <- function(object, ...) {
  if (!inherits(object, "cpaic_maic")) {
    stop("effective_sample_size() is only defined for cMAIC fits.",
         call. = FALSE)
  }
  object$ess
}
