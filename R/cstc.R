# Component STC: anchored simulated treatment comparison + component bridge ---

#' Anchored STC for one IPD study
#'
#' Fits a family-appropriate outcome regression with treatment main
#' effects, prognostic main effects, and treatment-by-effect-modifier
#' interactions, where the effect modifiers are centered at the target
#' population means. The treatment coefficient is then the
#' population-adjusted (anchored) contrast versus the reference arm in the
#' target population (the interaction terms vanish at centered = 0).
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
  # Center effect modifiers at the target population; only EM centering
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
                       function(z) sum(z != 0))
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
#' interactions; the effect modifiers are centered at a common `target`
#' population so the treatment coefficient is the anchored,
#' population-adjusted contrast in that population. These adjusted
#' contrasts replace the corresponding unadjusted aggregate contrasts and
#' [cnma_bridge()] combines them through the additive component model.
#'
#' Unlike [cmaic()] (reweighting) this is the regression-adjustment route.
#' The reported treatment coefficient is the *conditional* effect at the
#' target effect-modifier means (not a marginal standardization); for
#' collapsible measures the two coincide. It is implemented natively here
#' because `mlumr::stc()` targets the *unanchored* two-trial case; the link
#' and standard-error machinery is adapted from that package.
#'
#' @param network A [cpaic_network()] object that includes IPD.
#' @param target Named numeric vector (or list / one-row data frame) of
#'   target-population means for the effect modifiers.
#' @param effect_modifiers Covariates that interact with treatment
#'   (centered at `target`). Defaults to all IPD covariates.
#' @param prognostics Covariates included as main effects only. Defaults to
#'   the effect modifiers (so each enters as main effect + interaction).
#' @param common,random Passed to [cnma_bridge()].
#'
#' @section What the two-stage bridge does and does not adjust:
#' Only the edges carrying individual patient data are population-adjusted to the
#' target. Every aggregate-only edge keeps its published, study-population
#' contrast, and the additive bridge then combines all edges as if they estimated
#' the same component effects. Under effect modification they do not: an aggregate
#' edge estimates its contrast in *its own* trial population, while the adjusted
#' IPD edge estimates it at the target. The two agree only when the aggregate
#' populations resemble the target, or when the components on those edges are not
#' effect-modified. Treat a cross-network contrast that leans on aggregate-only
#' edges as adjusted for the IPD part alone, and prefer [cmlnmr()], which carries
#' the component by effect-modifier interactions through the whole network and so
#' adjusts every edge to the same target population coherently.
#'
#' @return An object of class `cpaic_stc` (and `cpaic_fit`).
#' @seealso [cmaic()], [cnma_bridge()]
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                      family = "binomial", ipd_covariates = "x1",
#'                      inactive = "Placebo")
#' fit <- cstc(net, target = c(x1 = 0), effect_modifiers = "x1")
#' relative_effects(fit)
#' additivity_test(fit)
#' @export
cstc <- function(network, target, effect_modifiers = NULL,
                 prognostics = NULL, common = FALSE, random = TRUE) {
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
  adj <- vector("list", length(info$studies))

  for (i in seq_along(info$studies)) {
    s <- info$studies[i]
    ipd_s <- network$ipd[as.character(network$ipd[[info$study]]) == s, ,
                         drop = FALSE]
    agd_s <- agd[as.character(agd[[cols$studlab]]) == s, , drop = FALSE]
    arms <- unique(as.character(ipd_s[[info$trt]]))
    if (length(arms) < 2L) {
      stop("IPD study '", s, "' has a single arm; cstc() needs a within-study ",
           "contrast (at least two arms). A one-arm study carries no anchored ",
           "treatment effect.", call. = FALSE)
    }
    if (length(arms) > 2L) {
      stop("IPD study '", s, "' has ", length(arms), " arms; cstc() ",
           "supports two-arm IPD studies in this version.", call. = FALSE)
    }
    ref_arm <- if (nrow(agd_s)) {
      as.character(agd_s[[cols$treat2]][1])
    } else if (network$reference %in% arms) {
      network$reference
    } else {
      sort(arms)[1]
    }
    if (!ref_arm %in% arms) ref_arm <- sort(arms)[1]

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
      network = network
    ),
    class = c("cpaic_stc", "cpaic_fit")
  )
}

#' @export
print.cpaic_stc <- function(x, ...) {
  cat("cpaic: component STC (anchored; IPD edges adjusted to target)\n")
  cat("  Diagnostic two-stage bridge: aggregate-only edges keep their own study\n",
      "  population, so the pooled result is only partially target-adjusted.\n",
      "  Prefer cmlnmr() for a coherent single-target synthesis.\n", sep = "")
  cat("  Effect modifiers (x treatment): ",
      paste(x$effect_modifiers, collapse = ", "), "\n", sep = "")
  cat("  Prognostic main effects:        ",
      paste(x$prognostics, collapse = ", "), "\n", sep = "")
  cat("\n")
  print(x$bridge)
  invisible(x)
}
