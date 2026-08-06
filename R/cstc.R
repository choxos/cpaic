# Component STC: anchored simulated treatment comparison + component bridge ---

#' Anchored STC for one IPD study
#'
#' Fits a family-appropriate outcome regression with treatment main
#' effects, prognostic main effects, and treatment-by-effect-modifier
#' interactions, where the effect modifiers are centered at target means. The
#' treatment coefficient is the anchored average conditional link-scale
#' contrast versus the reference arm at those means.
#' @noRd
.cpaic_stc_one_study <- function(ipd_s, info, family, ref_arm, target_mean,
                                 effect_modifiers, prognostics, sm,
                                 outcome_args, study_id) {
  arm_col <- info$trt
  out_col <- info$outcome
  d <- ipd_s

  # Use an internal, collision-proof factor name so a user covariate cannot
  # silently be treated as the treatment factor.
  if (".cpaic_arm" %in% names(d)) {
    stop("cstc(): the reserved column `.cpaic_arm` collides with an IPD column ",
         "in study '", study_id, "'; rename it.", call. = FALSE)
  }
  # Center effect modifiers at the target means; only EM centering
  # affects the treatment coefficient (no treatment x prognostic terms).
  for (em in effect_modifiers) d[[em]] <- d[[em]] - target_mean[[em]]
  d$.cpaic_arm <- stats::relevel(factor(d[[arm_col]]), ref = ref_arm)

  # Quote every raw column name so a non-syntactic covariate name cannot break
  # or, worse, silently alter the model formula.
  bq <- function(x) sprintf("`%s`", x)
  rhs <- ".cpaic_arm"
  if (length(prognostics)) {
    rhs <- paste(rhs, "+", paste(bq(prognostics), collapse = " + "))
  }
  if (length(effect_modifiers)) {
    rhs <- paste(rhs, "+",
                 paste(sprintf(".cpaic_arm:%s", bq(effect_modifiers)),
                       collapse = " + "))
  }

  # Capture (do not suppress) fit warnings so separation / non-convergence can
  # be classified by the validity gate.
  warn <- character(0)
  capture <- function(expr) withCallingHandlers(
    expr,
    warning = function(cnd) {
      warn <<- c(warn, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })

  n_events <- NULL
  if (family == "survival") {
    f <- stats::as.formula(paste0(
      "survival::Surv(", bq(outcome_args$time), ", ", bq(outcome_args$status),
      ") ~ ", rhs))
    fit <- capture(survival::coxph(f, data = d))
    n_events <- tapply(d[[outcome_args$status]], d$.cpaic_arm,
                       function(z) sum(z == 1L))
  } else {
    fam <- switch(family,
                  binomial = stats::binomial(),
                  gaussian = stats::gaussian(),
                  poisson  = stats::poisson(),
                  stop("unsupported family: ", family, call. = FALSE))
    if (family == "poisson" && !is.null(outcome_args$exposure)) {
      f <- stats::as.formula(paste0(bq(out_col), " ~ ", rhs, " + offset(log(",
                                    bq(outcome_args$exposure), "))"))
    } else {
      f <- stats::as.formula(paste0(bq(out_col), " ~ ", rhs))
    }
    fit <- capture(stats::glm(f, family = fam, data = d))
  }
  cf <- stats::coef(fit)
  V <- stats::vcov(fit)

  # The treatment main-effect coefficients are known exactly from the arm factor
  # levels; deriving them by regex risks matching an interaction term or a
  # covariate whose name happens to start with the internal prefix.
  arm_coef <- paste0(".cpaic_arm", setdiff(levels(d$.cpaic_arm), ref_arm))
  probs <- .cpaic_regression_problems(
    fit, family, expected_terms = arm_coef, n_events = n_events, warn = warn)
  if (length(probs)) .cpaic_stop_invalid_edge("cstc()", study_id, probs)

  est <- cf[arm_coef]
  se <- sqrt(diag(V)[arm_coef])
  arms_non_ref <- sub("^\\.cpaic_arm", "", names(est))

  list(
    contrasts = data.frame(
      treat1 = arms_non_ref, treat2 = ref_arm,
      TE = unname(est), seTE = unname(se),
      stringsAsFactors = FALSE),
    n = nrow(ipd_s), fit = fit
  )
}

#' Component simulated treatment comparison (cSTC)
#'
#' Anchored STC generalized to a (possibly disconnected) component network.
#' For each IPD study an outcome regression is fitted with treatment
#' main effects, prognostic main effects, and treatment-by-effect-modifier
#' interactions. The effect modifiers are centered at common target means, so
#' the treatment coefficient is the anchored average conditional link-scale
#' contrast at those means. These adjusted
#' contrasts replace the corresponding unadjusted aggregate contrasts and
#' [cnma_bridge()] combines them through the additive component model.
#'
#' Unlike [cmaic()] (reweighting) this is the regression-adjustment route.
#' The reported treatment coefficient is the *conditional* effect at the
#' target effect-modifier means. Equivalently, under the fitted linear
#' interaction model this is the average conditional link-scale effect at the
#' supplied target means, not a marginal standardization. It is implemented natively here
#' because the `stc()` function in the mlumr package targets the *unanchored*
#' two-trial case; the link and standard-error machinery is adapted from that
#' package. (Written without the double-colon form on purpose: mlumr is not a
#' dependency of cpaic, and the documentation site resolves a qualified
#' package-and-function reference by loading that package.)
#'
#' @param network A [cpaic_network()] object that includes IPD.
#' @param target Named numeric vector (or list / one-row data frame) of target
#'   means for the effect modifiers.
#' @param effect_modifiers Covariates that interact with treatment
#'   (centered at `target`). Defaults to all IPD covariates.
#' @param prognostics Covariates included as main effects only. Defaults to
#'   the effect modifiers (so each enters as main effect + interaction).
#' @param reference Optional anchor (comparator) arm to use in every IPD study
#'   in which it appears, instead of inferring it from the aggregate row order.
#' @param common,random Passed to [cnma_bridge()].
#' @param allow_experimental_bridge Logical. The default `FALSE` stops when
#'   aggregate-only edges would be combined with target-adjusted IPD edges.
#'   Set `TRUE` only for explicitly exploratory sensitivity work; the fit
#'   records the retained edges and the reason the bridge is approximate.
#' @param allow_ipd_only_studies Logical. The default `FALSE` requires every
#'   IPD study to match exactly one aggregate two-arm edge. Set `TRUE` to append
#'   an IPD-derived edge that has no aggregate row. Such additions are recorded
#'   in the returned fit.
#'
#' @section What the two-stage bridge does and does not adjust:
#' Only the edges carrying individual patient data are population-adjusted to the
#' target means. Every aggregate-only edge keeps its published study-specific
#' contrast, and the additive bridge then combines all edges as if they estimated
#' the same component effects. Under effect modification they do not: an aggregate
#' edge estimates its contrast in *its own* trial population, while the adjusted
#' IPD edge estimates it at the target. The two agree only when the aggregate
#' populations resemble the target, or when the components on those edges are not
#' effect-modified. Treat a cross-network contrast that leans on aggregate-only
#' edges as adjusted for the IPD part alone. Prefer [cmlnmr()] for a joint model
#' whose average conditional link-scale outputs are explicitly evaluated at
#' common target effect-modifier means.
#'
#' @return An object of class `cpaic_stc` (and `cpaic_fit`).
#' @seealso [cmaic()], [cnma_bridge()]
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                      family = "binomial", ipd_covariates = "x1",
#'                      inactive = "Placebo")
#' fit <- cstc(net, target = c(x1 = 0), effect_modifiers = "x1",
#'              allow_experimental_bridge = TRUE)
#' relative_effects(fit)
#' additivity_test(fit)
#' @export
cstc <- function(network, target, effect_modifiers = NULL,
                 prognostics = NULL, common = FALSE, random = TRUE,
                 reference = NULL, allow_experimental_bridge = FALSE,
                 allow_ipd_only_studies = FALSE) {
  stopifnot(inherits(network, "cpaic_network"))
  if (is.null(network$ipd)) {
    stop("`network` has no IPD; cstc() requires individual patient data.",
         call. = FALSE)
  }
  info <- network$ipd_info
  family <- network$family
  if (is.null(effect_modifiers)) effect_modifiers <- info$covariates
  if (is.null(prognostics)) prognostics <- effect_modifiers
  # Interaction hierarchy: an effect modifier that interacts with treatment must
  # also enter as a main effect, or the treatment-by-modifier interaction has no
  # corresponding main effect and the fit is not interpretable.
  prognostics <- union(prognostics, effect_modifiers)
  target <- as.list(target)
  target_mean <- target[effect_modifiers]
  if (any(vapply(target_mean, is.null, logical(1)))) {
    stop("`target` must supply a mean for every effect modifier: ",
         paste(effect_modifiers, collapse = ", "), call. = FALSE)
  }
  if (!all(vapply(target_mean, function(v)
           is.numeric(v) && length(v) == 1L && is.finite(v), logical(1)))) {
    stop("`target` values must be finite numeric scalars.", call. = FALSE)
  }

  outcome_args <- list(time = network$cols$ipd_time,
                       status = network$cols$ipd_status,
                       exposure = network$cols$ipd_exposure)
  agd <- network$agd
  cols <- network$cols
  plan <- .cpaic_two_stage_plan(
    network, reference, "cstc()",
    allow_ipd_only_studies = allow_ipd_only_studies)
  bridge_validity <- .cpaic_two_stage_bridge_gate(
    agd, plan$adjusted, cols, "cstc()", family,
    allow_experimental_bridge = allow_experimental_bridge)
  bridge_validity$ipd_only_studies <- plan$ipd_only_studies
  adj <- vector("list", length(plan$studies))

  for (i in seq_along(plan$studies)) {
    study_plan <- plan$studies[[i]]
    s <- study_plan$study
    ipd_s <- study_plan$ipd
    ref_arm <- study_plan$reference

    res <- .cpaic_stc_one_study(ipd_s, info, family, ref_arm, target_mean,
                                effect_modifiers, prognostics, network$sm,
                                outcome_args, study_id = s)
    res$contrasts[[cols$studlab]] <- s
    adj[[i]] <- res$contrasts
  }

  adj_df <- do.call(rbind, adj)
  agd2 <- .cpaic_replace_contrasts(agd, adj_df, cols)
  net2 <- network
  net2$agd <- agd2
  bridge <- cnma_bridge(net2, common = common, random = random)

  structure(
    list(
      bridge = bridge,
      components = bridge$components,
      target = target_mean,
      effect_modifiers = effect_modifiers,
      prognostics = prognostics,
      method = "cSTC",
      adjusted_contrasts = adj_df,
      ipd_only_studies = plan$ipd_only_studies,
      bridge_validity = bridge_validity,
      network = network
    ),
    class = c("cpaic_stc", "cpaic_fit")
  )
}

#' @export
print.cpaic_stc <- function(x, ...) {
  cat("cpaic: component STC (anchored; IPD edges at target means)\n")
  validity <- x$bridge_validity
  if (is.null(validity)) {
    cat("  Two-stage bridge gate: NOT RECORDED\n",
        "  Refit before interpreting this result for decision-grade use.\n",
        sep = "")
  } else if (isTRUE(validity$decision_grade)) {
    cat("  Two-stage bridge gate: PASSED (decision-grade eligibility only)\n",
        "  This gate result does not establish overall model validity.\n",
        sep = "")
  } else {
    override <- if (isTRUE(validity$experimental_override)) {
      "EXPERIMENTAL OVERRIDE ACTIVE"
    } else {
      "FAILED"
    }
    cat("  Two-stage bridge gate: ", override, "\n", sep = "")
    cat("  This bridge is not decision-grade; interpret it only as exploratory\n",
        "  sensitivity output.\n", sep = "")
    if (length(validity$reasons)) {
      cat("  Reasons:\n")
      for (reason in validity$reasons) cat("    - ", reason, "\n", sep = "")
    }
    retained <- validity$retained_aggregate_edges
    cols <- x$network$cols
    required <- unlist(cols[c("studlab", "treat1", "treat2")], use.names = FALSE)
    if (is.data.frame(retained) && nrow(retained) &&
        length(required) == 3L && all(required %in% names(retained))) {
      labels <- paste0(
        as.character(retained[[cols$studlab]]), ": ",
        as.character(retained[[cols$treat1]]), " vs ",
        as.character(retained[[cols$treat2]]))
      cat("  Retained aggregate-only edges: ",
          paste(labels, collapse = "; "), "\n", sep = "")
    }
  }
  if (length(validity$ipd_only_studies)) {
    cat("  Explicitly appended IPD-only studies: ",
        paste(validity$ipd_only_studies, collapse = ", "), "\n", sep = "")
  }
  cat("  Effect modifiers (x treatment): ",
      paste(x$effect_modifiers, collapse = ", "), "\n", sep = "")
  cat("  Prognostic main effects:        ",
      paste(x$prognostics, collapse = ", "), "\n", sep = "")
  cat("\n")
  print(x$bridge)
  invisible(x)
}
