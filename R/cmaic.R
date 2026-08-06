# Component MAIC: target-matched, component-bridged indirect comparison ---------

#' Weighted within-study contrast for one IPD study
#'
#' Fits the family-appropriate weighted outcome model and returns the
#' estimated contrast of each non-reference arm versus the reference arm,
#' on the link / log scale used by `sm`.
#' @noRd
.cpaic_weighted_fit <- function(data, family, arm_col, ref_arm,
                                outcome_col, weights,
                                time_col = NULL, status_col = NULL,
                                exposure_col = NULL) {
  # Note: do not reorder `data` here; `weights` is aligned to the
  # incoming row order; `relevel()` fixes the reference arm regardless.
  arm <- stats::relevel(factor(data[[arm_col]]), ref = ref_arm)
  w <- weights
  # Capture (do not suppress) warnings so the validity gate can classify
  # separation and non-convergence rather than proceeding blindly.
  warn <- character(0)
  capture <- function(expr) withCallingHandlers(
    expr, warning = function(cnd) {
      warn <<- c(warn, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  n_events <- NULL

  if (family == "survival") {
    if (is.null(time_col) || is.null(status_col)) {
      stop("survival family needs `time` and `status` columns.", call. = FALSE)
    }
    fit <- capture(survival::coxph(
      survival::Surv(data[[time_col]], data[[status_col]]) ~ arm,
      weights = w, robust = TRUE))
    n_events <- tapply(data[[status_col]], arm, function(z) sum(z == 1L))
  } else {
    fam <- switch(family,
                  binomial = stats::binomial(),
                  gaussian = stats::gaussian(),
                  poisson  = stats::poisson(),
                  stop("unsupported family: ", family, call. = FALSE))
    y <- data[[outcome_col]]
    if (family == "poisson" && !is.null(exposure_col)) {
      if (any(!is.finite(data[[exposure_col]]) | data[[exposure_col]] <= 0)) {
        stop("Poisson IPD exposure must be positive and finite.",
             call. = FALSE)
      }
      off <- log(data[[exposure_col]])
      fit <- capture(stats::glm(y ~ arm + offset(off), family = fam, weights = w))
    } else {
      fit <- capture(stats::glm(y ~ arm, family = fam, weights = w))
    }
  }
  cf <- stats::coef(fit)
  cf <- cf[grepl("^arm", names(cf))]
  if (!length(cf)) {
    stop("No non-reference arm coefficient was estimated for an IPD study; ",
         "check the treatment/arm coding.", call. = FALSE)
  }
  names(cf) <- sub("^arm", "", names(cf))
  list(cf = cf, fit = fit, warn = warn, n_events = n_events)
}

# Thin wrapper for the bootstrap loop, where only the coefficient is needed.
#' @noRd
.cpaic_weighted_contrast <- function(...) .cpaic_weighted_fit(...)$cf

#' Center effect modifiers on target means
#' @noRd
.cpaic_center <- function(data, target_mean, target_sd = NULL) {
  ems <- names(target_mean)
  for (em in ems) {
    data[[paste0(em, "_CENTERED")]] <- data[[em]] - target_mean[[em]]
    if (!is.null(target_sd) && !is.null(target_sd[[em]]) &&
        !is.na(target_sd[[em]])) {
      data[[paste0(em, "_sq_CENTERED")]] <-
        data[[em]]^2 - (target_mean[[em]]^2 + target_sd[[em]]^2)
    }
  }
  data
}

#' Classify bootstrap validity failures
#' @noRd
.cpaic_bootstrap_problem_codes <- function(stage, problems) {
  vapply(problems, function(problem) {
    if (stage == "weight_validation") {
      if (grepl("non-finite or negative weights", problem, fixed = TRUE)) {
        return("weights_nonfinite_or_negative")
      }
      if (grepl("weights sum to zero", problem, fixed = TRUE)) {
        return("weights_zero_sum")
      }
      if (grepl("non-finite or zero ESS", problem, fixed = TRUE)) {
        return("ess_nonfinite_or_zero")
      }
      if (grepl("optimizer did not converge", problem, fixed = TRUE)) {
        return("weight_optimizer_nonconvergence")
      }
      if (grepl("moment balance", problem, fixed = TRUE)) {
        return("moment_imbalance")
      }
      return("weight_validation_failure")
    }
    if (stage == "regression_validation") {
      if (grepl("non-finite coefficient", problem, fixed = TRUE)) {
        return("coefficient_nonfinite")
      }
      if (grepl("missing treatment coefficient", problem, fixed = TRUE)) {
        return("treatment_coefficient_missing")
      }
      if (grepl("no covariance matrix", problem, fixed = TRUE)) {
        return("covariance_missing")
      }
      if (grepl("degenerate covariance", problem, fixed = TRUE)) {
        return("treatment_covariance_degenerate")
      }
      if (grepl("did not converge", problem, fixed = TRUE)) {
        return("outcome_model_nonconvergence")
      }
      if (grepl("aliased", problem, fixed = TRUE)) {
        return("treatment_term_aliased")
      }
      if (grepl("zero events", problem, fixed = TRUE)) {
        return("arm_zero_events")
      }
      if (grepl("separated / non-identified", problem, fixed = TRUE)) {
        return("treatment_effect_nonidentified")
      }
      if (grepl("fit warning", problem, fixed = TRUE)) {
        return("outcome_model_warning")
      }
      return("regression_validation_failure")
    }
    paste0(stage, "_failure")
  }, character(1))
}

#' Build bootstrap failure records
#' @noRd
.cpaic_bootstrap_failure_rows <- function(replicate, stage, problems,
                                          codes = NULL) {
  problems <- as.character(problems)
  if (!length(problems)) problems <- "unspecified bootstrap failure"
  if (is.null(codes)) {
    codes <- .cpaic_bootstrap_problem_codes(stage, problems)
  }
  data.frame(
    replicate = rep.int(as.integer(replicate), length(problems)),
    stage = rep.int(stage, length(problems)),
    reason_code = as.character(codes),
    reason = problems,
    stringsAsFactors = FALSE
  )
}

#' Summarize cMAIC bootstrap validity
#' @noRd
.cpaic_bootstrap_validation <- function(boot, accepted, failure_records,
                                        arms_non_ref, ref_arm, n_boot,
                                        min_boot_success) {
  failures <- do.call(rbind, Filter(Negate(is.null), failure_records))
  if (is.null(failures)) {
    failures <- data.frame(
      replicate = integer(), stage = character(), reason_code = character(),
      reason = character(), stringsAsFactors = FALSE)
  }

  # Defensive accounting invariant: every rejected row must be recorded even
  # if a future branch forgets to register a specific reason.
  recorded <- unique(failures$replicate)
  unrecorded <- setdiff(which(!accepted), recorded)
  if (length(unrecorded)) {
    failure_records[unrecorded] <- lapply(unrecorded, function(replicate) {
      rbind(
        failure_records[[replicate]],
        .cpaic_bootstrap_failure_rows(
          replicate, "internal", "bootstrap replicate was not accepted",
          "unclassified_rejection")
      )
    })
    failures <- do.call(rbind, Filter(Negate(is.null), failure_records))
  }

  se <- apply(boot, 2, stats::sd, na.rm = TRUE)
  n_ok <- colSums(is.finite(boot))
  success_fraction <- n_ok / n_boot
  se_mcse <- se / sqrt(2 * pmax(n_ok - 1L, 1L))

  # At least 20 successful replicates gives a normal-theory relative Monte
  # Carlo error for the bootstrap SD of about 16 percent or less. When a caller
  # deliberately requests fewer than 20 replicates for a small test fixture, all
  # requested replicates must succeed. Production analyses should use far more
  # than this floor; the default remains 500.
  absolute_min_success <- min(20L, n_boot)
  proportion_min_success <- as.integer(ceiling(min_boot_success * n_boot))
  required_success <- max(absolute_min_success, proportion_min_success)

  bootstrap_summary <- data.frame(
    treat1 = arms_non_ref,
    treat2 = ref_arm,
    n_requested = rep.int(n_boot, length(arms_non_ref)),
    n_success = unname(n_ok[arms_non_ref]),
    success_fraction = unname(success_fraction[arms_non_ref]),
    required_success_count = rep.int(required_success,
                                     length(arms_non_ref)),
    required_success_fraction = rep.int(required_success / n_boot,
                                        length(arms_non_ref)),
    bootstrap_se = unname(se[arms_non_ref]),
    bootstrap_se_mcse_normal_approx = unname(se_mcse[arms_non_ref]),
    stringsAsFactors = FALSE
  )
  failure_table <- if (nrow(failures)) {
    tab <- stats::aggregate(
      replicate ~ stage + reason_code,
      data = unique(failures[c("replicate", "stage", "reason_code")]),
      FUN = length)
    names(tab)[names(tab) == "replicate"] <- "n_replicates"
    tab$fraction_of_requested <- tab$n_replicates / n_boot
    tab
  } else {
    data.frame(
      stage = character(), reason_code = character(),
      n_replicates = integer(), fraction_of_requested = numeric(),
      stringsAsFactors = FALSE)
  }

  list(
    se = se,
    n_ok = n_ok,
    required_success = required_success,
    failures = failures,
    bootstrap_summary = bootstrap_summary,
    bootstrap_failure_table = failure_table,
    bootstrap_mcse_method = paste0(
      "Normal-theory Monte Carlo standard error of the bootstrap SD: ",
      "bootstrap_se / sqrt(2 * (n_success - 1))."),
    bootstrap_success_rule = paste0(
      "Each contrast requires at least max(ceiling(min_boot_success * n_boot), ",
      "min(20, n_boot)) successful replicates. When n_boot is below 20, all ",
      "requested replicates must succeed.")
  )
}

#' Stop for insufficient cMAIC bootstrap validity
#' @noRd
.cpaic_stop_bootstrap_validation_failure <- function(study_id, validation,
                                                      n_boot,
                                                      min_boot_success, boot) {
  if (!any(validation$n_ok < validation$required_success)) {
    return(invisible(NULL))
  }

  message <- paste0(
    "cmaic(): the population-adjusted fit for study '", study_id,
    "' is not usable and would corrupt the component bridge:\n  - only ",
    min(validation$n_ok), " of ", n_boot, " bootstrap replicates succeeded ",
    "(required ", validation$required_success, ": max of ceiling(",
    min_boot_success, " * ", n_boot, ") and min(20, ", n_boot,
    ")); the standard error is unreliable (poor overlap or separation in ",
    "the resamples)\nFix or remove this study; an invalid edge is not ",
    "silently dropped.")
  stop(structure(
    list(
      message = message,
      call = NULL,
      study = study_id,
      bootstrap_draws = boot,
      bootstrap_summary = validation$bootstrap_summary,
      bootstrap_failures = validation$failures,
      bootstrap_failure_table = validation$bootstrap_failure_table,
      bootstrap_mcse_method = validation$bootstrap_mcse_method,
      bootstrap_success_rule = validation$bootstrap_success_rule
    ),
    class = c("cpaic_bootstrap_error", "error", "condition")))
}

#' MAIC for one IPD study: weights + adjusted contrast(s) with bootstrap SE
#' @noRd
.cpaic_maic_one_study <- function(ipd_s, info, family, ref_arm,
                                  target_mean, target_sd, em_centered_cols,
                                  n_boot, min_boot_success, outcome_args,
                                  study_id) {
  arm_col <- info$trt
  out_col <- info$outcome

  centered <- .cpaic_center(ipd_s, target_mean, target_sd)
  wfit <- suppressMessages(maicplus::estimate_weights(
    centered, centered_colnames = em_centered_cols, boot_strata = arm_col))
  w <- wfit$data$weights
  ess <- wfit$ess

  # Gate the weight solution: finite positive weights, usable ESS, optimizer
  # convergence, and (critically) that the weights actually balanced the
  # requested moments. An unbalanced or degenerate solution must not become an
  # adjusted edge.
  wp <- .cpaic_weight_problems(w, ess, centered, em_centered_cols,
                               opt = wfit$opt)
  if (length(wp)) .cpaic_stop_invalid_edge("cmaic()", study_id, wp)

  pf <- .cpaic_weighted_fit(
    centered, family, arm_col, ref_arm, out_col, weights = w,
    time_col = outcome_args$time, status_col = outcome_args$status,
    exposure_col = outcome_args$exposure)
  arm_terms <- paste0("arm", setdiff(unique(as.character(centered[[arm_col]])),
                                     ref_arm))
  rp <- .cpaic_regression_problems(pf$fit, family, expected_terms = arm_terms,
                                   n_events = pf$n_events, warn = pf$warn)
  if (length(rp)) .cpaic_stop_invalid_edge("cmaic()", study_id, rp)
  point <- pf$cf

  # Stratified (by arm) bootstrap: re-estimate weights and refit, propagating
  # both the weighting and the outcome-model uncertainty.
  arms_non_ref <- names(point)
  boot <- matrix(NA_real_, nrow = n_boot, ncol = length(arms_non_ref),
                 dimnames = list(as.character(seq_len(n_boot)), arms_non_ref))
  strata <- split(seq_len(nrow(centered)), centered[[arm_col]])
  accepted <- rep(FALSE, n_boot)
  failure_records <- vector("list", n_boot)

  record_failure <- function(replicate, stage, problems, codes = NULL) {
    failure_records[[replicate]] <<- rbind(
      failure_records[[replicate]],
      .cpaic_bootstrap_failure_rows(replicate, stage, problems, codes)
    )
    invisible(NULL)
  }

  for (b in seq_len(n_boot)) {
    idx <- unlist(lapply(strata, function(ii) sample(ii, length(ii),
                                                     replace = TRUE)),
                  use.names = FALSE)
    db <- centered[idx, , drop = FALSE]
    weight_error <- NULL
    wfit_b <- tryCatch(
      suppressMessages(maicplus::estimate_weights(
        db, centered_colnames = em_centered_cols,
        boot_strata = arm_col)),
      error = function(e) {
        weight_error <<- conditionMessage(e)
        NULL
      })
    if (is.null(wfit_b)) {
      record_failure(b, "weight_estimation", weight_error,
                     "weight_estimation_error")
      next
    }
    # Hold a resampled weight solution to the SAME validity gate as the point
    # estimate: convergence, positivity, a usable ESS, and the moment balance
    # actually achieved. A failing replicate is a failed replicate.
    weight_gate <- tryCatch({
      wb <- wfit_b$data$weights
      list(
        ok = TRUE,
        weights = wb,
        problems = .cpaic_weight_problems(
          wb, wfit_b$ess, db, em_centered_cols, opt = wfit_b$opt)
      )
    }, error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    })
    if (!weight_gate$ok) {
      record_failure(b, "weight_validation", weight_gate$error,
                     "weight_validation_error")
      next
    }
    if (length(weight_gate$problems)) {
      record_failure(b, "weight_validation", weight_gate$problems)
      next
    }
    wb <- weight_gate$weights
    outcome_error <- NULL
    pf_b <- tryCatch(
      .cpaic_weighted_fit(db, family, arm_col, ref_arm, out_col,
                          weights = wb, time_col = outcome_args$time,
                          status_col = outcome_args$status,
                          exposure_col = outcome_args$exposure),
      error = function(e) {
        outcome_error <<- conditionMessage(e)
        NULL
      })
    if (is.null(pf_b)) {
      record_failure(b, "outcome_fit", outcome_error, "outcome_fit_error")
      next
    }
    regression_gate <- tryCatch({
      list(
        ok = TRUE,
        problems = .cpaic_regression_problems(
          pf_b$fit, family, expected_terms = arm_terms,
          n_events = pf_b$n_events, warn = pf_b$warn)
      )
    }, error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    })
    if (!regression_gate$ok) {
      record_failure(b, "regression_validation", regression_gate$error,
                     "regression_validation_error")
      next
    }
    if (length(regression_gate$problems)) {
      record_failure(b, "regression_validation", regression_gate$problems)
      next
    }
    cb <- pf_b$cf
    if (!setequal(names(cb), arms_non_ref) || any(!is.finite(cb))) {
      record_failure(
        b, "coefficient_validation",
        "validated outcome fit did not return one finite coefficient per contrast",
        "coefficient_set_invalid")
      next
    }
    boot[b, arms_non_ref] <- cb[arms_non_ref]
    accepted[b] <- TRUE
  }

  validation <- .cpaic_bootstrap_validation(
    boot, accepted, failure_records, arms_non_ref, ref_arm, n_boot,
    min_boot_success)
  .cpaic_stop_bootstrap_validation_failure(
    study_id, validation, n_boot, min_boot_success, boot)

  list(
    contrasts = data.frame(
      treat1 = arms_non_ref, treat2 = ref_arm,
      TE = unname(point), seTE = unname(validation$se[arms_non_ref]),
      stringsAsFactors = FALSE),
    ess = ess, weights = w, n = nrow(ipd_s),
    diagnostics = .cpaic_weight_diagnostics(w, centered, em_centered_cols),
    bootstrap_draws = boot,
    bootstrap_summary = validation$bootstrap_summary,
    bootstrap_failures = validation$failures,
    bootstrap_failure_table = validation$bootstrap_failure_table,
    bootstrap_mcse_method = validation$bootstrap_mcse_method,
    bootstrap_success_rule = validation$bootstrap_success_rule
  )
}

#' Component matching-adjusted indirect comparison (cMAIC)
#'
#' Anchored MAIC generalized to a (possibly disconnected) component
#' network. Each IPD study is reweighted with [maicplus::estimate_weights()]
#' so that its requested effect-modifier moments match a common `target`;
#' the resulting target-matched within-study contrasts
#' (with bootstrap standard errors that propagate the weighting
#' uncertainty) then replace the corresponding unadjusted aggregate
#' contrasts. Finally [cnma_bridge()] combines all contrasts through the
#' additive component model. The bridge is gated because retained aggregate
#' edges and nonlinear marginal effects can make that synthesis incoherent.
#'
#' @section What the two-stage bridge does and does not adjust:
#' Only the edges carrying individual patient data are population-adjusted to the
#' target moments. Every aggregate-only edge keeps its published study-specific
#' contrast, and the additive bridge then combines all edges as if they estimated
#' the same component effects. Under effect modification they do not: an aggregate
#' edge estimates its contrast in *its own* trial population, while the reweighted
#' IPD edge estimates it at the target. The two agree only when the aggregate
#' populations resemble the target, or when the components on those edges are not
#' effect-modified. Treat a cross-network contrast that leans on aggregate-only
#' edges as adjusted for the IPD part alone. Prefer [cmlnmr()] for a joint model
#' whose average conditional link-scale outputs are explicitly evaluated at
#' common target effect-modifier means.
#'
#' @section Non-collapsibility and the additive model:
#' cMAIC returns a **marginal** effect in the reweighted IPD sample, and the additive
#' component model assumes effects add. On a non-collapsible scale (the odds
#' ratio, the hazard ratio) **marginal effects do not add**, even when every
#' conditional effect does. In one simulated target population the marginal
#' log-odds ratios satisfied
#' `marginal(A) + marginal(B) = 0.6615` while `marginal(A+B) = 0.6411`; the
#' additive model is simply false on that scale. cMAIC therefore carries an
#' **irreducible approximation error** that survives perfect matching and
#' infinite sample size. Its size is problem-specific and cannot be assumed
#' negligible.
#'
#' Marginal component effects are not *generally* additive; they add exactly when
#' the standardized treatment effects remain affine in the component design.
#' Additivity is therefore a property of the conditional link scale that the
#' marginal scale inherits only approximately, and the error does not vanish with
#' sample size. Where it is material, [cstc()] or [cmlnmr()], which target a
#' conditional effect and inherit additivity exactly, are preferable. Note also
#' that the two-stage route combines a conditional adjusted edge with aggregate
#' edges reported on a marginal scale, so it should be regarded as approximate.
#'
#' @param network A [cpaic_network()] object that includes IPD.
#' @param target Named numeric vector (or one-row data frame / list) giving
#'   target means of the effect modifiers.
#' @param effect_modifiers Character vector of covariates to match on
#'   (defaults to all IPD covariates). Matching only on effect modifiers is
#'   the anchored-MAIC convention.
#' @param target_sd Optional named numeric vector of target standard
#'   deviations; when supplied, second moments are matched as well.
#' @param n_boot Number of bootstrap resamples for the adjusted-contrast
#'   standard errors. Default `500`.
#' @param min_boot_success Minimum fraction of bootstrap resamples that must
#'   succeed for a contrast. The enforced count is
#'   `max(ceiling(min_boot_success * n_boot), min(20, n_boot))`, so a run with
#'   fewer than 20 requested resamples requires every resample to succeed.
#'   Below this threshold the edge is rejected rather than given a fragile
#'   standard error from a selected subset. Default `0.8`.
#' @param reference Optional anchor (comparator) arm to use in every IPD study
#'   in which it appears, instead of inferring it from the aggregate row order.
#' @param seed Optional RNG seed for reproducible bootstrap. The caller's global
#'   RNG state is restored on exit, so calling `cmaic()` does not perturb a
#'   downstream random stream.
#' @param common,random Passed to [cnma_bridge()].
#' @param allow_experimental_bridge Logical. The default `FALSE` stops when
#'   aggregate-only edges would be combined with target-matched IPD edges, or
#'   when a non-Gaussian cMAIC contrast would be forced through an additive
#'   component model. Set `TRUE` only for explicitly exploratory sensitivity
#'   work; the fit records the exact approximation reasons.
#' @param allow_ipd_only_studies Logical. The default `FALSE` requires every
#'   IPD study to match exactly one aggregate two-arm edge. Set `TRUE` to append
#'   an IPD-derived edge that has no aggregate row. Such additions are recorded
#'   in the returned fit.
#'
#' @return An object of class `cpaic_maic` (also inheriting `cpaic_bridge`
#'   structure via `$bridge`), with the bridged fit, per-study effective
#'   sample sizes, and the target moments. Bootstrap diagnostic fields include
#'   `$bootstrap_draws`, `$bootstrap_summary`, `$bootstrap_failures`,
#'   `$bootstrap_failure_table`, `$bootstrap_mcse_method`, and
#'   `$bootstrap_success_rule`. A threshold failure raises a
#'   `cpaic_bootstrap_error` condition carrying the same diagnostic information.
#' @seealso [cstc()], [cnma_bridge()]
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                      family = "binomial", ipd_covariates = "x1",
#'                      inactive = "Placebo")
#' \donttest{
#' fit <- cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
#'              n_boot = 100, seed = 1,
#'              allow_experimental_bridge = TRUE)
#' relative_effects(fit)
#' effective_sample_size(fit)
#' }
#' @export
cmaic <- function(network, target, effect_modifiers = NULL, target_sd = NULL,
                  n_boot = 500, min_boot_success = 0.8, seed = NULL,
                  common = FALSE, random = TRUE, reference = NULL,
                  allow_experimental_bridge = FALSE,
                  allow_ipd_only_studies = FALSE) {
  stopifnot(inherits(network, "cpaic_network"))
  if (is.null(network$ipd)) {
    stop("`network` has no IPD; cmaic() requires individual patient data.",
         call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot < 2L || n_boot != as.integer(n_boot)) {
    stop("`n_boot` must be an integer >= 2.", call. = FALSE)
  }
  n_boot <- as.integer(n_boot)
  if (!is.numeric(min_boot_success) || length(min_boot_success) != 1L ||
      !is.finite(min_boot_success) || min_boot_success <= 0 ||
      min_boot_success > 1) {
    stop("`min_boot_success` must be a fraction in (0, 1].", call. = FALSE)
  }
  # Local RNG scope: restore the caller's stream on exit so chaining fits or
  # downstream Monte Carlo is not perturbed by our bootstrap resampling.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = .GlobalEnv))
        rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  info <- network$ipd_info
  family <- network$family
  if (is.null(effect_modifiers)) effect_modifiers <- info$covariates
  target <- as.list(target)
  target_mean <- target[effect_modifiers]
  if (anyNA(names(target_mean)) || any(vapply(target_mean, is.null, logical(1)))) {
    stop("`target` must supply a mean for every effect modifier: ",
         paste(effect_modifiers, collapse = ", "), call. = FALSE)
  }
  if (!all(vapply(target_mean, function(v)
           is.numeric(v) && length(v) == 1L && is.finite(v), logical(1)))) {
    stop("`target` values must be finite numeric scalars.", call. = FALSE)
  }
  em_centered_cols <- paste0(effect_modifiers, "_CENTERED")
  # Match second moments too when target_sd is supplied (otherwise the
  # _sq_CENTERED columns built by .cpaic_center() are never matched).
  if (!is.null(target_sd)) {
    target_sd <- as.list(target_sd)
    if (any(vapply(target_sd, function(v)
            !is.null(v) && (!is.numeric(v) || length(v) != 1L ||
                            !is.finite(v) || v < 0), logical(1)))) {
      stop("`target_sd` values must be finite and non-negative.",
           call. = FALSE)
    }
    has_sd <- vapply(effect_modifiers, function(e)
      !is.null(target_sd[[e]]) && is.finite(target_sd[[e]]), logical(1))
    em_centered_cols <- c(em_centered_cols,
                          paste0(effect_modifiers[has_sd], "_sq_CENTERED"))
  }

  outcome_args <- list(time = network$cols$ipd_time,
                       status = network$cols$ipd_status,
                       exposure = network$cols$ipd_exposure)

  agd <- network$agd
  cols <- network$cols
  plan <- .cpaic_two_stage_plan(
    network, reference, "cmaic()",
    allow_ipd_only_studies = allow_ipd_only_studies)
  bridge_validity <- .cpaic_two_stage_bridge_gate(
    agd, plan$adjusted, cols, "cmaic()", family,
    allow_experimental_bridge = allow_experimental_bridge)
  bridge_validity$ipd_only_studies <- plan$ipd_only_studies
  adj <- vector("list", length(plan$studies))
  wdiag <- vector("list", length(plan$studies))
  bootstrap_draws <- setNames(vector("list", length(plan$studies)),
                              vapply(plan$studies, `[[`, character(1), "study"))
  bootstrap_summary <- vector("list", length(plan$studies))
  bootstrap_failures <- vector("list", length(plan$studies))
  bootstrap_failure_table <- vector("list", length(plan$studies))
  ess <- setNames(numeric(length(info$studies)), info$studies)

  for (i in seq_along(plan$studies)) {
    study_plan <- plan$studies[[i]]
    s <- study_plan$study
    ipd_s <- study_plan$ipd
    ref_arm <- study_plan$reference

    res <- .cpaic_maic_one_study(
      ipd_s, info, family, ref_arm, target_mean, target_sd,
      em_centered_cols, n_boot, min_boot_success, outcome_args, study_id = s)
    res$contrasts[[cols$studlab]] <- s
    adj[[i]] <- res$contrasts
    ess[s] <- res$ess
    res$diagnostics <- cbind(study = s, res$diagnostics)
    wdiag[[i]] <- res$diagnostics
    bootstrap_draws[[i]] <- res$bootstrap_draws
    bootstrap_summary[[i]] <- cbind(
      study = s, res$bootstrap_summary, stringsAsFactors = FALSE)
    bootstrap_failures[[i]] <- cbind(
      study = rep.int(s, nrow(res$bootstrap_failures)),
      res$bootstrap_failures, stringsAsFactors = FALSE)
    bootstrap_failure_table[[i]] <- cbind(
      study = rep.int(s, nrow(res$bootstrap_failure_table)),
      res$bootstrap_failure_table, stringsAsFactors = FALSE)
  }

  adj_df <- do.call(rbind, adj)
  wdiag_df <- do.call(rbind, wdiag)
  bootstrap_summary_df <- do.call(rbind, bootstrap_summary)
  bootstrap_failures_df <- do.call(rbind, bootstrap_failures)
  bootstrap_failure_table_df <- do.call(rbind, bootstrap_failure_table)
  rownames(bootstrap_summary_df) <- NULL
  rownames(bootstrap_failures_df) <- NULL
  rownames(bootstrap_failure_table_df) <- NULL
  agd2 <- .cpaic_replace_contrasts(agd, adj_df, cols)

  net2 <- network
  net2$agd <- agd2
  bridge <- cnma_bridge(net2, common = common, random = random)

  structure(
    list(
      bridge = bridge,
      components = bridge$components,
      ess = ess,
      weight_diagnostics = wdiag_df,
      target = target_mean,
      effect_modifiers = effect_modifiers,
      n_boot = n_boot,
      bootstrap_draws = bootstrap_draws,
      bootstrap_summary = bootstrap_summary_df,
      bootstrap_failures = bootstrap_failures_df,
      bootstrap_failure_table = bootstrap_failure_table_df,
      bootstrap_mcse_method = res$bootstrap_mcse_method,
      bootstrap_success_rule = res$bootstrap_success_rule,
      method = "cMAIC",
      adjusted_contrasts = adj_df,
      ipd_only_studies = plan$ipd_only_studies,
      bridge_validity = bridge_validity,
      network = network
    ),
    class = c("cpaic_maic", "cpaic_fit")
  )
}

#' Weight-quality diagnostics for a cMAIC fit
#'
#' Per IPD study, the effective sample size, weight-entropy efficiency,
#' coefficient of variation, largest normalized weight, mass in the top 5% of
#' weights, and the largest residual effect-modifier imbalance after weighting.
#' A high maximum weight or low entropy efficiency signals a few dominant
#' individuals, which the effective sample size alone can hide.
#'
#' @param object A [cmaic()] fit.
#' @return A data frame, one row per IPD study.
#' @seealso [cmaic()], [effective_sample_size()]
#' @export
weight_diagnostics <- function(object) {
  if (!inherits(object, "cpaic_maic") || is.null(object$weight_diagnostics)) {
    stop("`object` must be a cmaic() fit.", call. = FALSE)
  }
  object$weight_diagnostics
}

#' Set one cell to a value, preserving the column's type (expanding factor
#' levels when needed) so an appended row does not silently coerce a column.
#' @noRd
.cpaic_set_cell <- function(col, value) {
  if (is.factor(col)) {
    factor(value, levels = union(levels(col), as.character(value)))
  } else if (is.numeric(col)) {
    suppressWarnings(as.numeric(value))
  } else {
    value
  }
}

#' rbind one edge onto `agd`, harmonizing factor levels first.
#' @noRd
.cpaic_rbind_edge <- function(agd, newrow) {
  for (nm in names(agd)) {
    if (is.factor(agd[[nm]]) && is.factor(newrow[[nm]])) {
      lv <- union(levels(agd[[nm]]), levels(newrow[[nm]]))
      agd[[nm]] <- factor(as.character(agd[[nm]]), levels = lv)
      newrow[[nm]] <- factor(as.character(newrow[[nm]]), levels = lv)
    }
  }
  rbind(agd, newrow)
}

#' Replace (or append) aggregate contrasts with population-adjusted ones
#'
#' Matching is orientation-insensitive (a study-`{treat1,treat2}` pair keyed
#' on the unordered treatment set), so an aggregate row recorded in the
#' opposite direction is *replaced* (with the sign flipped) rather than
#' appended, which would double-count the study. Both frames are required to
#' hold a unique key: two aggregate rows for the same study and pair would each
#' be overwritten with the same adjusted estimate and counted twice. An appended
#' edge is built from a typed-NA prototype, never cloned from an arbitrary
#' existing row, so it carries no unrelated metadata.
#' @noRd
.cpaic_replace_contrasts <- function(agd, adj_df, cols) {
  ukey <- function(sl, a, b) {
    paste(sl, pmin(a, b), pmax(a, b), sep = "\r")
  }
  agd_key <- ukey(as.character(agd[[cols$studlab]]),
                  as.character(agd[[cols$treat1]]),
                  as.character(agd[[cols$treat2]]))
  if (anyDuplicated(agd_key)) {
    dup <- unique(agd_key[duplicated(agd_key)])
    stop("The aggregate data contain ", length(dup), " duplicate ",
         "{study, treatment-pair} row(s); a duplicated edge would be counted ",
         "twice by the additive bridge. De-duplicate `agd` before adjusting.",
         call. = FALSE)
  }
  adj_key <- ukey(as.character(adj_df[[cols$studlab]]),
                  as.character(adj_df$treat1), as.character(adj_df$treat2))
  if (anyDuplicated(adj_key)) {
    stop("Duplicate adjusted contrasts for the same {study, treatment-pair}.",
         call. = FALSE)
  }

  # One typed-NA prototype row, so an appended contrast never inherits unrelated
  # metadata (covariate summaries, sample sizes, custom columns) from row 1.
  proto <- agd[1, , drop = FALSE]
  proto[] <- lapply(proto, function(col) col[NA_integer_])

  for (j in seq_len(nrow(adj_df))) {
    sl <- as.character(adj_df[[cols$studlab]][j])
    t1 <- as.character(adj_df$treat1[j])
    t2 <- as.character(adj_df$treat2[j])
    hit <- which(agd_key == ukey(sl, t1, t2))
    if (length(hit) == 1L) {
      h <- hit
      same_dir <- as.character(agd[[cols$treat1]][h]) == t1 &&
        as.character(agd[[cols$treat2]][h]) == t2
      agd[[cols$TE]][h]   <- if (same_dir) adj_df$TE[j] else -adj_df$TE[j]
      agd[[cols$seTE]][h] <- adj_df$seTE[j]
    } else {
      newrow <- proto
      newrow[[cols$studlab]] <- .cpaic_set_cell(newrow[[cols$studlab]], sl)
      newrow[[cols$treat1]]  <- .cpaic_set_cell(newrow[[cols$treat1]], t1)
      newrow[[cols$treat2]]  <- .cpaic_set_cell(newrow[[cols$treat2]], t2)
      newrow[[cols$TE]]      <- adj_df$TE[j]
      newrow[[cols$seTE]]    <- adj_df$seTE[j]
      agd <- .cpaic_rbind_edge(agd, newrow)
    }
  }
  agd
}

#' @export
component_effects.cpaic_fit <- function(object, newdata = NULL, ...) {
  object$components
}

#' @export
print.cpaic_maic <- function(x, ...) {
  cat("cpaic: component MAIC (anchored; IPD edges matched to target moments)\n")
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
  cat("  Effect modifiers matched: ",
      paste(x$effect_modifiers, collapse = ", "), "\n", sep = "")
  cat("  Effective sample sizes (per IPD study):\n")
  for (s in names(x$ess)) {
    cat("    ", s, ": ESS = ", round(x$ess[s], 1), "\n", sep = "")
  }
  cat("\n")
  print(x$bridge)
  invisible(x)
}
