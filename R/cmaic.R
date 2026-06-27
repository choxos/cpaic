# Component MAIC: population-adjusted, component-bridged indirect comparison ---

#' Weighted within-study contrast for one IPD study
#'
#' Fits the family-appropriate weighted outcome model and returns the
#' estimated contrast of each non-reference arm versus the reference arm,
#' on the link / log scale used by `sm`.
#' @noRd
.cpaic_weighted_contrast <- function(data, family, arm_col, ref_arm,
                                     outcome_col, weights,
                                     time_col = NULL, status_col = NULL,
                                     exposure_col = NULL) {
  # Note: do not reorder `data` here -- `weights` is aligned to the
  # incoming row order; `relevel()` fixes the reference arm regardless.
  arm <- stats::relevel(factor(data[[arm_col]]), ref = ref_arm)
  w <- weights

  if (family == "survival") {
    if (is.null(time_col) || is.null(status_col)) {
      stop("survival family needs `time` and `status` columns.", call. = FALSE)
    }
    fit <- survival::coxph(
      survival::Surv(data[[time_col]], data[[status_col]]) ~ arm,
      weights = w, robust = TRUE
    )
    cf <- stats::coef(fit)
  } else {
    fam <- switch(family,
                  binomial = stats::binomial(),
                  gaussian = stats::gaussian(),
                  poisson  = stats::poisson(),
                  stop("unsupported family: ", family, call. = FALSE))
    y <- data[[outcome_col]]
    if (family == "poisson" && !is.null(exposure_col)) {
      off <- log(data[[exposure_col]])
      fit <- stats::glm(y ~ arm + offset(off), family = fam, weights = w)
    } else {
      fit <- suppressWarnings(stats::glm(y ~ arm, family = fam, weights = w))
    }
    cf <- stats::coef(fit)
  }
  cf <- cf[grepl("^arm", names(cf))]
  names(cf) <- sub("^arm", "", names(cf))
  cf
}

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
.cpaic_maic_one_study <- function(ipd_s, info, family, sm, ref_arm,
                                  target_mean, target_sd, em_centered_cols,
                                  n_boot, outcome_args) {
  arm_col <- info$trt
  out_col <- info$outcome

  centered <- .cpaic_center(ipd_s, target_mean, target_sd)
  wfit <- suppressMessages(maicplus::estimate_weights(
    centered, centered_colnames = em_centered_cols, boot_strata = arm_col))
  w <- wfit$data$weights
  ess <- wfit$ess

  point <- .cpaic_weighted_contrast(
    centered, family, arm_col, ref_arm, out_col, weights = w,
    time_col = outcome_args$time, status_col = outcome_args$status,
    exposure_col = outcome_args$exposure)

  # Stratified (by arm) bootstrap: re-estimate weights and refit, propagating
  # both the weighting and the outcome-model uncertainty.
  arms_non_ref <- names(point)
  boot <- matrix(NA_real_, nrow = n_boot, ncol = length(arms_non_ref),
                 dimnames = list(NULL, arms_non_ref))
  strata <- split(seq_len(nrow(centered)), centered[[arm_col]])
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
      error = function(e) rep(NA_real_, length(arms_non_ref)))
    boot[b, names(cb)] <- cb
  }
  se <- apply(boot, 2, stats::sd, na.rm = TRUE)

  list(
    contrasts = data.frame(
      treat1 = arms_non_ref, treat2 = ref_arm,
      TE = unname(point), seTE = unname(se[arms_non_ref]),
      stringsAsFactors = FALSE),
    ess = ess, weights = w, n = nrow(ipd_s)
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
#' @param seed Optional RNG seed for reproducible bootstrap.
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
                  n_boot = 500, seed = NULL, common = FALSE, random = TRUE) {
  stopifnot(inherits(network, "cpaic_network"))
  if (is.null(network$ipd)) {
    stop("`network` has no IPD; cmaic() requires individual patient data.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
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
  ess <- setNames(numeric(length(info$studies)), info$studies)

  for (i in seq_along(info$studies)) {
    s <- info$studies[i]
    ipd_s <- network$ipd[as.character(network$ipd[[info$study]]) == s, ,
                         drop = FALSE]
    # Reference (anchor) arm: the comparator of this study's AgD row if
    # present, else the network reference if it is an arm, else first arm.
    agd_s <- agd[as.character(agd[[cols$studlab]]) == s, , drop = FALSE]
    arms <- unique(as.character(ipd_s[[info$trt]]))
    if (length(arms) > 2L) {
      stop("IPD study '", s, "' has ", length(arms), " arms; cmaic() ",
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

    res <- .cpaic_maic_one_study(
      ipd_s, info, family, network$sm, ref_arm, target_mean, target_sd,
      em_centered_cols, n_boot, outcome_args)
    res$contrasts[[cols$studlab]] <- s
    adj[[i]] <- res$contrasts
    ess[s] <- res$ess
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
      ess = ess,
      target = target_mean,
      effect_modifiers = effect_modifiers,
      n_boot = n_boot,
      method = "cMAIC",
      network = network
    ),
    class = c("cpaic_maic", "cpaic_fit")
  )
}

#' Replace (or append) aggregate contrasts with population-adjusted ones
#' @noRd
.cpaic_replace_contrasts <- function(agd, adj_df, cols) {
  key <- function(df, t1, t2, sl) paste(df[[sl]], df[[t1]], df[[t2]], sep = "\r")
  agd_key <- key(agd, cols$treat1, cols$treat2, cols$studlab)
  adj_key <- paste(adj_df[[cols$studlab]], adj_df$treat1, adj_df$treat2,
                   sep = "\r")
  for (j in seq_len(nrow(adj_df))) {
    hit <- which(agd_key == adj_key[j])
    if (length(hit)) {
      agd[[cols$TE]][hit]   <- adj_df$TE[j]
      agd[[cols$seTE]][hit] <- adj_df$seTE[j]
    } else {
      newrow <- agd[1, , drop = FALSE]
      newrow[[cols$studlab]] <- adj_df[[cols$studlab]][j]
      newrow[[cols$treat1]]  <- adj_df$treat1[j]
      newrow[[cols$treat2]]  <- adj_df$treat2[j]
      newrow[[cols$TE]]      <- adj_df$TE[j]
      newrow[[cols$seTE]]    <- adj_df$seTE[j]
      agd <- rbind(agd, newrow)
    }
  }
  agd
}

#' @export
component_effects.cpaic_fit <- function(object, ...) object$components

#' @export
print.cpaic_maic <- function(x, ...) {
  cat("cpaic: component MAIC (anchored, population-adjusted)\n")
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
