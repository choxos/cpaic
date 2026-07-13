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
       "(e.g. LOO).", call. = FALSE)
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
    cat("  Additivity restrictions (Q.diff): not available -- no standard NMA\n",
        "    is estimable on a disconnected network.\n", sep = "")
  }
  cat("  Note: neither statistic tests whether component effects are constant\n",
      "  ACROSS sub-networks, which is the assumption that bridges the gap.\n",
      "  That assumption is untestable from the data and must be justified\n",
      "  clinically.\n", sep = "")
  invisible(x)
}

#' Effective sample sizes from a cMAIC fit
#'
#' @param object A `cpaic_maic` object.
#' @param ... Unused.
#' @return A named numeric vector of effective sample sizes per IPD study.
#' @export
effective_sample_size <- function(object, ...) {
  if (!inherits(object, "cpaic_maic")) {
    stop("effective_sample_size() is only defined for cMAIC fits.",
         call. = FALSE)
  }
  object$ess
}
