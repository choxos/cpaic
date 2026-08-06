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

#' Weight-quality diagnostics for one MAIC study
#'
#' All computed from the achieved weights, so they are free. ESS alone hides a
#' few extreme weights or residual imbalance; these expose them.
#' @noRd
.cpaic_weight_diagnostics <- function(w, centered, centered_cols) {
  n <- length(w)
  sw <- sum(w)
  p <- w / sw
  ord <- sort(w, decreasing = TRUE)
  top <- sum(ord[seq_len(max(1L, ceiling(0.05 * n)))]) / sw
  present <- intersect(centered_cols, names(centered))
  bal <- if (length(present)) {
    max(abs(vapply(present, function(cc) sum(w * centered[[cc]]) / sw,
                   numeric(1))))
  } else NA_real_
  pnz <- p[p > 0]
  data.frame(
    ess = sw^2 / sum(w^2),
    n = n,
    entropy_eff = (-sum(pnz * log(pnz))) / log(n),
    cv = stats::sd(w) / mean(w),
    max_weight = max(w) / sw,
    top5pct_mass = top,
    max_abs_balance = bal,
    row.names = NULL, stringsAsFactors = FALSE)
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

#' Validate that a two-stage IPD fit replaces a complete two-arm study edge
#'
#' A study with more aggregate contrasts than the two IPD arms is a partially
#' observed multi-arm study. Replacing one contrast while retaining the others
#' discards their joint covariance and creates a hybrid study that neither the
#' outcome model nor the aggregate network model represents.
#' @noRd
.cpaic_validate_two_stage_edge <- function(agd_s, cols, arms, study_id,
                                           method) {
  if (nrow(agd_s) > 1L) {
    stop(method, ": study '", study_id, "' has partial IPD representation of ",
         "a multi-arm aggregate study. All within-study contrasts share ",
         "sampling covariance, so one pair cannot be replaced while the ",
         "remaining aggregate pairs are retained. Supply complete IPD for a ",
         "supported joint multi-arm model, or use a separate two-arm study.",
         call. = FALSE)
  }
  if (nrow(agd_s) == 1L) {
    pair <- as.character(agd_s[1, c(cols$treat1, cols$treat2)])
    if (!setequal(pair, arms)) {
      stop(method, ": study '", study_id, "' has IPD arms {",
           paste(sort(arms), collapse = ", "), "}, but its aggregate edge is {",
           paste(sort(pair), collapse = ", "), "}. The adjusted contrast cannot ",
           "replace a different treatment pair.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Preflight and plan each two-stage IPD contrast before any model is fitted
#' @noRd
.cpaic_two_stage_plan <- function(network, reference, method,
                                  allow_ipd_only_studies = FALSE) {
  if (!is.logical(allow_ipd_only_studies) ||
      length(allow_ipd_only_studies) != 1L || is.na(allow_ipd_only_studies)) {
    stop("`allow_ipd_only_studies` must be TRUE or FALSE.", call. = FALSE)
  }
  info <- network$ipd_info
  agd <- network$agd
  cols <- network$cols
  studies <- lapply(info$studies, function(s) {
    ipd_s <- network$ipd[
      as.character(network$ipd[[info$study]]) == s, , drop = FALSE]
    agd_s <- agd[as.character(agd[[cols$studlab]]) == s, , drop = FALSE]
    arms <- unique(as.character(ipd_s[[info$trt]]))
    if (length(arms) < 2L) {
      stop("IPD study '", s, "' has a single arm; ", method,
           " needs a within-study contrast (at least two arms).", call. = FALSE)
    }
    if (length(arms) > 2L) {
      stop("IPD study '", s, "' has ", length(arms), " arms; ", method,
           " supports two-arm IPD studies in this version.", call. = FALSE)
    }
    if (!nrow(agd_s) && !allow_ipd_only_studies) {
      stop(method, ": IPD study '", s, "' has no aggregate row to replace. ",
           "Appending an IPD-only edge changes the evidence set. Set ",
           "`allow_ipd_only_studies = TRUE` explicitly to permit and record ",
           "that addition.", call. = FALSE)
    }
    .cpaic_validate_two_stage_edge(agd_s, cols, arms, s, method)
    ref_arm <- if (!is.null(reference) && reference %in% arms) {
      reference
    } else if (nrow(agd_s)) {
      as.character(agd_s[[cols$treat2]][1])
    } else if (network$reference %in% arms) {
      network$reference
    } else {
      sort(arms)[1]
    }
    if (!ref_arm %in% arms) ref_arm <- sort(arms)[1]
    list(study = s, ipd = ipd_s, arms = arms, reference = ref_arm,
         ipd_only = !nrow(agd_s))
  })

  planned <- do.call(rbind, lapply(studies, function(x) {
    out <- data.frame(
      treat1 = setdiff(x$arms, x$reference), treat2 = x$reference,
      stringsAsFactors = FALSE)
    out[[cols$studlab]] <- x$study
    out
  }))
  rownames(planned) <- NULL
  ipd_only <- vapply(
    studies, function(x) if (isTRUE(x$ipd_only)) x$study else NA_character_,
    character(1))
  ipd_only <- unname(ipd_only[!is.na(ipd_only)])
  list(
    studies = studies,
    adjusted = planned,
    ipd_only_studies = ipd_only
  )
}

#' Stable unordered key for a study-specific treatment contrast
#' @noRd
.cpaic_edge_key <- function(study, treat1, treat2) {
  paste(as.character(study), pmin(as.character(treat1), as.character(treat2)),
        pmax(as.character(treat1), as.character(treat2)), sep = "\r")
}

#' Gate methodologically approximate two-stage component bridges
#'
#' The two-stage estimators adjust only the IPD-bearing contrasts. Any retained
#' aggregate contrast remains tied to its own study population. In addition,
#' cMAIC produces marginal contrasts, which are not generally additive on a
#' nonlinear link scale. Both approximations require explicit opt-in.
#' @noRd
.cpaic_two_stage_bridge_gate <- function(agd, adjusted, cols, method, family,
                                         allow_experimental_bridge) {
  if (!is.logical(allow_experimental_bridge) ||
      length(allow_experimental_bridge) != 1L ||
      is.na(allow_experimental_bridge)) {
    stop("`allow_experimental_bridge` must be TRUE or FALSE.", call. = FALSE)
  }

  agd_key <- .cpaic_edge_key(agd[[cols$studlab]], agd[[cols$treat1]],
                             agd[[cols$treat2]])
  adjusted_key <- .cpaic_edge_key(
    adjusted[[cols$studlab]], adjusted$treat1, adjusted$treat2)
  retained <- agd[!agd_key %in% adjusted_key, , drop = FALSE]

  reasons <- character(0)
  if (nrow(retained)) {
    labels <- paste0(
      as.character(retained[[cols$studlab]]), ": ",
      as.character(retained[[cols$treat1]]), " vs ",
      as.character(retained[[cols$treat2]]))
    reasons <- c(reasons, paste0(
      "retained aggregate-only edge(s) remain in their own study populations: ",
      paste(labels, collapse = "; ")))
  }
  if (identical(method, "cmaic()") && !identical(family, "gaussian")) {
    reasons <- c(reasons, paste0(
      "cMAIC estimates marginal ", family, " contrasts, which are not generally ",
      "additive in the component design on a nonlinear link scale"))
  }

  if (length(reasons)) {
    msg <- paste0(
      method, " cannot form a decision-grade component bridge:\n  - ",
      paste(reasons, collapse = "\n  - "),
      "\nUse cmlnmr() for a joint model, restrict the analysis to a design in which ",
      "every edge is adjusted and the estimand is additive, or set ",
      "`allow_experimental_bridge = TRUE` only for explicitly exploratory ",
      "sensitivity work.")
    if (!allow_experimental_bridge) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }

  list(
    decision_grade = !length(reasons),
    experimental_override = length(reasons) && allow_experimental_bridge,
    reasons = reasons,
    retained_aggregate_edges = retained
  )
}
