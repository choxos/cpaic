# Component MAIC: population-adjusted, component-bridged indirect comparison ---

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
    n_events <- tapply(data[[status_col]], arm, function(z) sum(z != 0))
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

#' Center effect modifiers on the target population
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
                 dimnames = list(NULL, arms_non_ref))
  strata <- split(seq_len(nrow(centered)), centered[[arm_col]])
  logfam <- family %in% c("binomial", "poisson", "survival")
  for (b in seq_len(n_boot)) {
    idx <- unlist(lapply(strata, function(ii) sample(ii, length(ii),
                                                     replace = TRUE)),
                  use.names = FALSE)
    db <- centered[idx, , drop = FALSE]
    wb <- tryCatch(
      suppressMessages(maicplus::estimate_weights(
        db, centered_colnames = em_centered_cols,
        boot_strata = arm_col))$data$weights,
      error = function(e) NULL)
    if (is.null(wb)) next
    cb <- tryCatch(
      .cpaic_weighted_contrast(db, family, arm_col, ref_arm, out_col,
                               weights = wb, time_col = outcome_args$time,
                               status_col = outcome_args$status,
                               exposure_col = outcome_args$exposure),
      error = function(e) stats::setNames(rep(NA_real_, length(arms_non_ref)),
                                          arms_non_ref))
    # A separated / degenerate resample yields a finite but absurd coefficient;
    # treat it as a failed replicate so it neither inflates the bootstrap SE nor
    # counts toward the success threshold. On a bounded (log / hazard) link this
    # is flagged only when the resample is BOTH huge and far from the validated
    # point estimate, so a genuinely large but well-identified effect is kept.
    pt <- point[names(cb)]
    cb[!is.finite(cb) |
       (logfam & abs(cb) > 30 & abs(cb - pt) > 20)] <- NA_real_
    boot[b, names(cb)] <- cb
  }
  se <- apply(boot, 2, stats::sd, na.rm = TRUE)
  n_ok <- colSums(!is.na(boot))
  # Fail closed: too few successful replicates means the SE cannot be trusted,
  # which usually signals poor overlap or separation in the resamples. Abstain
  # rather than emit a fragile SE from a selected subset.
  if (any(n_ok < min_boot_success * n_boot)) {
    .cpaic_stop_invalid_edge("cmaic()", study_id, paste0(
      "only ", min(n_ok), " of ", n_boot, " bootstrap replicates succeeded ",
      "(threshold ", round(min_boot_success * n_boot), "); the standard error ",
      "is unreliable (poor overlap or separation in the resamples)"))
  }

  list(
    contrasts = data.frame(
      treat1 = arms_non_ref, treat2 = ref_arm,
      TE = unname(point), seTE = unname(se[arms_non_ref]),
      stringsAsFactors = FALSE),
    ess = ess, weights = w, n = nrow(ipd_s),
    diagnostics = .cpaic_weight_diagnostics(w, centered, em_centered_cols)
  )
}

#' Component matching-adjusted indirect comparison (cMAIC)
#'
#' Anchored MAIC generalized to a (possibly disconnected) component
#' network. Each IPD study is reweighted with [maicplus::estimate_weights()]
#' so that its effect-modifier distribution matches a common `target`
#' population; the resulting population-adjusted within-study contrasts
#' (with bootstrap standard errors that propagate the weighting
#' uncertainty) then replace the corresponding unadjusted aggregate
#' contrasts. Finally [cnma_bridge()] combines all contrasts through the
#' additive component model, yielding relative effects that are both
#' connected across sub-networks and adjusted to the target population.
#'
#' @section What the two-stage bridge does and does not adjust:
#' Only the edges carrying individual patient data are population-adjusted to the
#' target. Every aggregate-only edge keeps its published, study-population
#' contrast, and the additive bridge then combines all edges as if they estimated
#' the same component effects. Under effect modification they do not: an aggregate
#' edge estimates its contrast in *its own* trial population, while the reweighted
#' IPD edge estimates it at the target. The two agree only when the aggregate
#' populations resemble the target, or when the components on those edges are not
#' effect-modified. Treat a cross-network contrast that leans on aggregate-only
#' edges as adjusted for the IPD part alone, and prefer [cmlnmr()], which carries
#' the component by effect-modifier interactions through the whole network and so
#' adjusts every edge to the same target population coherently.
#'
#' @section Non-collapsibility and the additive model:
#' cMAIC returns a **marginal** effect in the target population, and the additive
#' component model assumes effects add. On a non-collapsible scale (the odds
#' ratio, the hazard ratio) **marginal effects do not add**, even when every
#' conditional effect does. In one simulated target population the marginal
#' log-odds ratios satisfied
#' `marginal(A) + marginal(B) = 0.6615` while `marginal(A+B) = 0.6411`; the
#' additive model is simply false on that scale. cMAIC therefore carries a small
#' **irreducible bias** (about +0.02 log-OR there) that survives perfect matching
#' and infinite sample size. It is small relative to a typical standard error
#' (about 0.25) but it does not vanish with more data.
#'
#' Marginal component effects are not *generally* additive; they add exactly when
#' the standardized treatment effects remain affine in the component design.
#' Additivity is therefore a property of the conditional link scale that the
#' marginal scale inherits only approximately, and the error does not vanish with
#' sample size. Where it is material, [cstc()] or [cmlnmr()], which target a
#' conditional effect and inherit additivity exactly, are preferable. Note also
#' that the two-stage route combines a conditional adjusted edge with aggregate
#' edges reported on a marginal scale, so it should be regarded as approximate. See
#' `documentation/validation/VALIDATION.md`.
#'
#' @param network A [cpaic_network()] object that includes IPD.
#' @param target Named numeric vector (or one-row data frame / list) giving
#'   the target-population means of the effect modifiers.
#' @param effect_modifiers Character vector of covariates to match on
#'   (defaults to all IPD covariates). Matching only on effect modifiers is
#'   the anchored-MAIC convention.
#' @param target_sd Optional named numeric vector of target standard
#'   deviations; when supplied, second moments are matched as well.
#' @param n_boot Number of bootstrap resamples for the adjusted-contrast
#'   standard errors. Default `500`.
#' @param min_boot_success Minimum fraction of bootstrap resamples that must
#'   succeed for a contrast; below this threshold the edge is rejected rather
#'   than given a fragile standard error from a selected subset. Default `0.8`.
#' @param reference Optional anchor (comparator) arm to use in every IPD study
#'   in which it appears, instead of inferring it from the aggregate row order.
#' @param seed Optional RNG seed for reproducible bootstrap. The caller's global
#'   RNG state is restored on exit, so calling `cmaic()` does not perturb a
#'   downstream random stream.
#' @param common,random Passed to [cnma_bridge()].
#'
#' @return An object of class `cpaic_maic` (also inheriting `cpaic_bridge`
#'   structure via `$bridge`), with the bridged fit, per-study effective
#'   sample sizes, and the target population.
#' @seealso [cstc()], [cnma_bridge()]
#' @examples
#' net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
#'                      family = "binomial", ipd_covariates = "x1",
#'                      inactive = "Placebo")
#' \donttest{
#' fit <- cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
#'              n_boot = 100, seed = 1)
#' relative_effects(fit)
#' effective_sample_size(fit)
#' }
#' @export
cmaic <- function(network, target, effect_modifiers = NULL, target_sd = NULL,
                  n_boot = 500, min_boot_success = 0.8, seed = NULL,
                  common = FALSE, random = TRUE, reference = NULL) {
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
  adj <- vector("list", length(info$studies))
  wdiag <- vector("list", length(info$studies))
  ess <- setNames(numeric(length(info$studies)), info$studies)

  for (i in seq_along(info$studies)) {
    s <- info$studies[i]
    ipd_s <- network$ipd[as.character(network$ipd[[info$study]]) == s, ,
                         drop = FALSE]
    # Reference (anchor) arm: the comparator of this study's AgD row if
    # present, else the network reference if it is an arm, else first arm.
    agd_s <- agd[as.character(agd[[cols$studlab]]) == s, , drop = FALSE]
    arms <- unique(as.character(ipd_s[[info$trt]]))
    if (length(arms) < 2L) {
      stop("IPD study '", s, "' has a single arm; cmaic() needs a within-study ",
           "contrast (at least two arms).", call. = FALSE)
    }
    if (length(arms) > 2L) {
      stop("IPD study '", s, "' has ", length(arms), " arms; cmaic() ",
           "supports two-arm IPD studies in this version.", call. = FALSE)
    }
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

    res <- .cpaic_maic_one_study(
      ipd_s, info, family, ref_arm, target_mean, target_sd,
      em_centered_cols, n_boot, min_boot_success, outcome_args, study_id = s)
    res$contrasts[[cols$studlab]] <- s
    adj[[i]] <- res$contrasts
    ess[s] <- res$ess
    res$diagnostics <- cbind(study = s, res$diagnostics)
    wdiag[[i]] <- res$diagnostics
  }

  adj_df <- do.call(rbind, adj)
  wdiag_df <- do.call(rbind, wdiag)
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
      method = "cMAIC",
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
  cat("cpaic: component MAIC (anchored; IPD edges adjusted to target)\n")
  cat("  Diagnostic two-stage bridge: aggregate-only edges keep their own study\n",
      "  population, so the pooled result is only partially target-adjusted.\n",
      "  Prefer cmlnmr() for a coherent single-target synthesis.\n", sep = "")
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
