# Diagnostics: additivity, effective sample size ------------------------------

#' Test the additivity assumption of the component model
#'
#' Returns the Cochran Q statistic for the additive component network
#' meta-analysis (from [netmeta::discomb()]). A small p-value indicates
#' lack of fit of the additive model, i.e. evidence of component
#' interactions; consider adding interaction terms.
#'
#' @param object A `cpaic_bridge` or `cpaic_fit` object.
#' @param ... Unused.
#' @return A one-row data frame with `Q`, `df`, and `pval`.
#' @export
additivity_test <- function(object, ...) {
  UseMethod("additivity_test")
}

#' @export
additivity_test.cpaic_fit <- function(object, ...) {
  additivity_test(object$bridge, ...)
}

#' @export
additivity_test.cpaic_bridge <- function(object, ...) {
  out <- data.frame(Q = object$Q$Q, df = object$Q$df, pval = object$Q$pval,
                    row.names = NULL)
  class(out) <- c("cpaic_additivity", "data.frame")
  out
}

#' @export
print.cpaic_additivity <- function(x, ...) {
  cat("Additivity (Cochran Q) test of the component model\n")
  cat("  Q = ", round(x$Q, 3), ", df = ", x$df, ", p = ",
      format.pval(x$pval, digits = 3), "\n", sep = "")
  if (!is.na(x$pval) && x$pval < 0.05) {
    cat("  -> evidence against additivity; consider interaction terms.\n")
  }
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
