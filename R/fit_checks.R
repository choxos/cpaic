# Fail-closed validity gates for two-stage regression and weight fits ---------
#
# The two-stage methods (cstc, cmaic) turn a per-study regression or weighted
# regression into a TE/seTE that is inserted into the additive component bridge.
# An invalid fit (non-convergence, separation, rank deficiency, a degenerate
# covariance) that is passed through unchecked injects a meaningless contrast
# into `discomb()` and silently corrupts every cross-gap effect that leans on
# it. These helpers detect such fits so the caller can stop rather than proceed.
#
# `fit$converged` alone is not enough: complete or quasi-complete separation in
# a binomial or Poisson GLM routinely reports `converged = TRUE` with a finite
# but absurd coefficient and an enormous standard error. The checks below add
# separation, boundary-fitted-value, rank/aliasing, finite-covariance, and
# implausible-magnitude tests on top of convergence.

#' Problems with a fitted GLM / Cox regression that feeds an adjusted edge
#'
#' Returns a character vector of problems; an empty vector means the fit is
#' usable. `expected_terms` are the treatment-coefficient names that must be
#' present and finite. `n_events` (survival only) guards zero-event arms.
#' `warn` carries any warnings captured while fitting.
#' @noRd
.cpaic_regression_problems <- function(fit, family, expected_terms,
                                       n_events = NULL, warn = character(0)) {
  probs <- character(0)
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  logfam <- family %in% c("binomial", "poisson", "survival")

  if (is.null(cf) || any(!is.finite(cf))) {
    probs <- c(probs, "non-finite coefficient(s)")
  }
  # A planned treatment coefficient absent from the fit means the design dropped
  # it (aliasing / rank deficiency).
  if (!is.null(cf)) {
    miss <- setdiff(expected_terms, names(cf))
    if (length(miss)) {
      probs <- c(probs, paste0("missing treatment coefficient(s): ",
                               paste(miss, collapse = ", ")))
    }
  }
  # Covariance must be finite and positive FOR THE TREATMENT TERMS only. A
  # degenerate nuisance/prognostic direction must not by itself fail an
  # otherwise-usable treatment contrast.
  present <- if (is.null(V)) character(0) else intersect(expected_terms,
                                                         rownames(V))
  if (is.null(V)) {
    probs <- c(probs, "no covariance matrix")
  } else if (length(present)) {
    var_tr <- diag(V)[present]
    if (any(!is.finite(var_tr)) || any(var_tr <= 0)) {
      probs <- c(probs, "degenerate covariance for a treatment coefficient")
    }
  }

  if (family != "survival") {
    if (isFALSE(fit$converged)) probs <- c(probs, "GLM did not converge")
    aliased <- tryCatch(summary(fit)$aliased, error = function(e) NULL)
    if (!is.null(aliased)) {
      tr_alias <- intersect(names(aliased), expected_terms)
      if (length(tr_alias) && any(aliased[tr_alias])) {
        probs <- c(probs, "aliased (rank-deficient) treatment term")
      }
    }
  } else if (!is.null(n_events) && any(n_events == 0L)) {
    probs <- c(probs, "an arm has zero events")
  }

  # Separation / non-identification on a bounded (log / hazard) link: the Wald
  # SE of a treatment effect explodes (a finite but absurd estimate paired with
  # an enormous SE). This is NOT applied to gaussian, where a mean-difference SE
  # can be large on an unbounded outcome scale without indicating separation,
  # and it is judged per observation only through the treatment coefficients, so
  # a strong prognostic covariate pushing individual fitted values to a boundary
  # does not trip it.
  if (logfam && !is.null(cf) && length(present)) {
    se_tr <- sqrt(pmax(diag(V)[present], 0))
    cf_tr <- cf[present]
    if (any(is.finite(se_tr) & se_tr > 30) ||
        any(is.finite(cf_tr) & is.finite(se_tr) & abs(cf_tr) > 30 & se_tr > 5)) {
      probs <- c(probs,
                 "separated / non-identified treatment effect (degenerate SE)")
    }
  }

  if (length(warn)) {
    bad <- grep("converge|did not|infinite|singular", warn,
                ignore.case = TRUE, value = TRUE)
    if (length(bad)) probs <- c(probs, paste0("fit warning: ", unique(bad)))
  }
  unique(probs)
}

#' Problems with a MAIC weight solution that feeds an adjusted edge
#'
#' Checks finiteness, positivity, effective sample size, optional optimizer
#' convergence, and (most importantly) whether the weights actually achieved
#' the requested moment balance. `centered` holds the `_CENTERED` columns whose
#' weighted means must be near zero for the match to be valid.
#' @noRd
.cpaic_weight_problems <- function(w, ess, centered, centered_cols,
                                   opt = NULL, tol = 1e-3) {
  probs <- character(0)
  if (any(!is.finite(w)) || any(w < 0)) {
    probs <- c(probs, "non-finite or negative weights")
  }
  if (!isTRUE(sum(w) > 0)) probs <- c(probs, "weights sum to zero")
  if (!is.finite(ess) || ess <= 0) probs <- c(probs, "non-finite or zero ESS")
  if (!is.null(opt) && !is.null(opt$convergence) &&
      !isTRUE(opt$convergence == 0)) {
    probs <- c(probs, "weight optimizer did not converge")
  }
  if (!length(probs) && sum(w) > 0) {
    present <- intersect(centered_cols, names(centered))
    if (length(present)) {
      wm <- vapply(present, function(cc)
        sum(w * centered[[cc]]) / sum(w), numeric(1))
      if (any(abs(wm) > tol)) {
        probs <- c(probs, paste0(
          "weights did not achieve moment balance (max |weighted mean| = ",
          signif(max(abs(wm)), 3), " > ", tol, ")"))
      }
    }
  }
  unique(probs)
}

#' Stop with a structured message when an adjusted edge is invalid
#' @noRd
.cpaic_stop_invalid_edge <- function(method, study_id, problems) {
  stop(method, ": the population-adjusted fit for study '", study_id,
       "' is not usable and would corrupt the component bridge:\n  - ",
       paste(problems, collapse = "\n  - "),
       "\nFix or remove this study; an invalid edge is not silently dropped.",
       call. = FALSE)
}
