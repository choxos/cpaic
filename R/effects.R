# Post-fit relative effects and league tables ---------------------------------

#' Validate a reference treatment and confidence/credible level
#' @noRd
.cpaic_check_ref_level <- function(reference, trts, level) {
  if (!is.character(reference) || length(reference) != 1L ||
      !reference %in% trts) {
    stop("`reference` must be one of the network treatments: ",
         paste(trts, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(NULL)
}

#' Is the summary measure on a (natural) log scale?
#' @noRd
.is_log_sm <- function(sm) {
  toupper(sm) %in% c("OR", "RR", "HR", "IRR", "RRR", "ROR", "HRR")
}

#' Validate an explicitly named supported estimand
#' @noRd
.cpaic_match_estimand <- function(estimand, supported, where) {
  if (!is.character(estimand) || length(estimand) != 1L || is.na(estimand) ||
      !nzchar(estimand)) {
    stop("`estimand` must be the single supported value \"", supported,
         "\" for ", where, ".", call. = FALSE)
  }
  if (!identical(estimand, supported)) {
    extra <- if (identical(estimand, "marginal")) {
      " cpaic does not yet implement marginal standardization for this output."
    } else {
      ""
    }
    stop(where, " supports estimand = \"", supported, "\" only.", extra,
         call. = FALSE)
  }
  estimand
}

#' Resolve target effect-modifier means and enforce their support
#'
#' The current summary path evaluates a conditional link-scale contrast at the
#' supplied covariate means. Because that contrast is linear in the effect
#' modifiers, it equals the average conditional link-scale contrast. It is not
#' a marginal standardized effect for nonlinear links.
#' @noRd
.cpaic_target_x <- function(newdata, effect_modifiers, margins = NULL) {
  Q <- length(effect_modifiers)
  if (!Q) return(numeric(0))
  if (is.null(newdata)) {
    stop("Specify `newdata`: a one-row data frame giving the target ",
         "effect-modifier means (",
         paste(effect_modifiers, collapse = ", "), "). Relative effects from ",
         "a model with effect modifiers depend on these means, so there is ",
         "no target-free answer. Use, for example, ",
         "newdata = data.frame(", effect_modifiers[1], " = 0) to obtain the ",
         "effect at the covariate origin.", call. = FALSE)
  }
  newdata <- as.data.frame(newdata)
  if (nrow(newdata) != 1L) {
    stop("`newdata` must have exactly one row of target means.",
         call. = FALSE)
  }
  miss <- setdiff(effect_modifiers, names(newdata))
  if (length(miss)) {
    stop("`newdata` is missing effect modifier(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  bad_type <- effect_modifiers[!vapply(effect_modifiers, function(v) {
    is.numeric(newdata[[v]]) && !is.factor(newdata[[v]])
  }, logical(1))]
  if (length(bad_type)) {
    stop("`newdata` effect-modifier values must be numeric, not factor or ",
         "character: ", paste(bad_type, collapse = ", "), ".", call. = FALSE)
  }
  x <- vapply(effect_modifiers, function(v) as.double(newdata[[v]][1]),
              numeric(1))
  if (any(!is.finite(x))) {
    stop("`newdata` effect-modifier values must be finite.", call. = FALSE)
  }
  if (!is.null(margins)) {
    if (is.null(names(margins))) {
      if (length(margins) != Q) {
        stop("Internal error: unnamed `margins` must have one entry per ",
             "effect modifier.", call. = FALSE)
      }
      names(margins) <- effect_modifiers
    }
    margins <- margins[effect_modifiers]
    if (anyNA(margins)) {
      stop("Internal error: `margins` is missing an effect modifier.",
           call. = FALSE)
    }
    for (v in effect_modifiers) {
      value <- x[[v]]
      margin <- margins[[v]]
      if (identical(margin, "bernoulli") && (value < 0 || value > 1)) {
        stop("Target mean `", v, "` must lie in [0, 1] for a Bernoulli ",
             "margin. Fractional values are valid prevalences.", call. = FALSE)
      }
      if (margin %in% c("gamma", "lognormal") && value <= 0) {
        stop("Target mean `", v, "` must be positive for a ", margin,
             " margin.", call. = FALSE)
      }
      if (identical(margin, "beta") && (value <= 0 || value >= 1)) {
        stop("Target mean `", v, "` must lie in (0, 1) for a beta margin.",
             call. = FALSE)
      }
    }
  }
  x
}

#' Posterior draws of component effects at target effect-modifier means
#'
#' In the component-additive ML-NMR the treatment effect of treatment `t` in a
#' model evaluated at effect-modifier means `x` is
#' `theta_t(x) = C_t' (beta + gamma x)`, so `beta + gamma x` are the component
#' average conditional link-scale effects at those means. Reporting `beta`
#' alone gives the effect at the covariate origin.
#' @noRd
.cpaic_beta_at <- function(object, x) {
  B <- .cpaic_draws_matrix(object$fit, "beta")
  nC <- ncol(object$C.matrix)
  if (!length(x)) return(B)
  G <- .cpaic_draws_matrix(object$fit, "gamma")
  out <- B
  for (q in seq_along(x)) {
    cols <- paste0("gamma[", seq_len(nC), ",", q, "]")
    miss <- setdiff(cols, colnames(G))
    if (length(miss)) {
      stop("The fit does not contain the component x effect-modifier draws ",
           "needed to evaluate the target means.", call. = FALSE)
    }
    out <- out + G[, cols, drop = FALSE] * x[q]
  }
  out
}

#' Relative treatment effects from a cpaic fit
#'
#' Tidies the relative effects of the fitted model: every treatment versus a
#' chosen reference, or all pairwise comparisons. Effects are reported on the
#' natural scale of the summary measure (e.g. odds ratios) unless
#' `backtransf = FALSE`.
#'
#' Relative effects that the component design cannot uniquely identify (their
#' contrast vector lies outside the row space of `X = B C`) are returned as
#' `NA` rather than as pseudoinverse or prior-driven artefacts. See
#' [estimable_effects()].
#'
#' For [cmlnmr()] fits the model contains component x effect-modifier
#' interactions, so relative effects depend on the supplied covariate means:
#' `theta_t(x) = C_t' (beta + gamma x)`. Because this expression is linear in
#' `x`, evaluating it at `E[X]` gives the average conditional link-scale effect.
#' It is not a marginal standardized OR, RR, or HR.
#'
#' Two things this returns are narrower than they may look, both documented in
#' full under [cmlnmr()]. The value is the **average conditional link-scale**
#' contrast at target means `x`, not the marginal effect standardized over a population
#' with a distribution of covariates; on a non-collapsible scale (odds ratio,
#' hazard ratio) those differ, so `newdata = <a study's covariate means>` does not
#' give that study's population-average effect. And where the interactions are
#' informed only by aggregate arms, `x` is being applied to an ecological
#' gradient rather than to within-study effect modification; see
#' [estimable_effects_at()], whose `identified_by` column separates the two.
#'
#' @param object A fitted cpaic object (`cpaic_bridge`, `cpaic_maic`,
#'   `cpaic_stc`, or `cpaic_mlnmr`).
#' @param reference Reference treatment. Defaults to the network reference.
#' @param all_contrasts If `TRUE`, return all pairwise comparisons instead
#'   of versus the reference.
#' @param backtransf If `TRUE` (default) back-transform log-scale measures
#'   (OR/RR/HR/...) by exponentiating.
#' @param level Confidence level for the intervals. Default `0.95`.
#' @param newdata For [cmlnmr()] fits: a one-row data frame giving target
#'   effect-modifier means. Required when the model has effect modifiers.
#' @param estimand For [cmlnmr()] fits, the only implemented value is
#'   `"average_conditional_link"`. Requests for `"marginal"` fail explicitly.
#' @param ... Unused.
#'
#' @return A data frame with columns `treatment`, `comparator`, `estimate`,
#'   `estimate_link`, `se_link`, `lower`, `upper`, `scale`, and `z`/`p` for
#'   frequentist fits. `estimate_link` and `se_link` are always on the model's
#'   link scale. `estimate`, `lower`, and `upper` use the scale named in `scale`.
#'   For [cmlnmr()] (Bayesian) fits the intervals are credible intervals and
#'   the final column is `pr_gt0`, the posterior probability that the
#'   effect (on the link scale) exceeds zero, instead of `z`/`p`.
#' @export
relative_effects <- function(object, reference = NULL, all_contrasts = FALSE,
                             backtransf = TRUE, level = 0.95, newdata = NULL,
                             estimand = NULL, ...) {
  UseMethod("relative_effects")
}

#' @export
relative_effects.cpaic_fit <- function(object, ...) {
  out <- relative_effects(object$bridge, ...)
  attr(out, "bridge_validity") <- object$bridge_validity %||% NULL
  out
}

#' @export
relative_effects.cpaic_mlnmr <- function(object, reference = NULL,
                                         all_contrasts = FALSE,
                                         backtransf = TRUE, level = 0.95,
                                         newdata = NULL,
                                         estimand = "average_conditional_link",
                                         ...) {
  .cpaic_match_estimand(estimand, "average_conditional_link",
                        "relative_effects() for cmlnmr fits")
  C <- object$C.matrix
  trts <- rownames(C)
  if (is.null(reference)) reference <- object$reference
  .cpaic_check_ref_level(reference, trts, level)

  x <- .cpaic_target_x(newdata, object$effect_modifiers, object$margins)
  Beff <- .cpaic_beta_at(object, x)      # draws x components at target means
  Theta <- Beff %*% t(C)                 # draws x treatments
  colnames(Theta) <- trts

  a <- (1 - level) / 2
  logsm <- .is_log_sm(object$sm)
  output_scale <- if (backtransf && logsm) "natural" else "link"
  # Estimability of an average conditional contrast depends on the target
  # means: it needs the relevant rows of Gamma, not just of beta. Test the
  # augmented contrast (1, x) %x% m against the joint information design. The
  # basis records how strong the identification is: "exact" only for an
  # injective-link IPD contrast, "first-order screen" when the row-space check
  # passes but leans on aggregate arms or a survival baseline (and so can be
  # optimistic), "not identified" otherwise.
  Njoint <- .cpaic_null_space(object$joint_design)
  Nipd <- .cpaic_null_space(object$joint_design_ipd)
  ipd_can_be_exact <- !identical(object$family, "survival")
  basis_of <- function(v) {
    if (!.cpaic_in_rowspace(matrix(v, nrow = 1L), Njoint)) return("not identified")
    if (ipd_can_be_exact && .cpaic_in_rowspace(matrix(v, nrow = 1L), Nipd)) {
      return("exact")
    }
    "first-order screen"
  }

  build <- function(t1, t2) {
    v <- .cpaic_target_vec(C[t1, ] - C[t2, ], x)
    bas <- basis_of(v)
    if (bas == "not identified") {
      return(data.frame(treatment = t1, comparator = t2, estimate = NA_real_,
                        estimate_link = NA_real_, se_link = NA_real_,
                        lower = NA_real_, upper = NA_real_,
                        scale = output_scale, pr_gt0 = NA_real_, basis = bas,
                        stringsAsFactors = FALSE))
    }
    d <- Theta[, t1] - Theta[, t2]                 # link scale
    e <- if (backtransf && logsm) exp(d) else d     # reporting scale
    q <- stats::quantile(e, c(a, 1 - a), names = FALSE)
    data.frame(treatment = t1, comparator = t2, estimate = mean(e),
               estimate_link = mean(d), se_link = stats::sd(d),
               lower = q[1], upper = q[2], scale = output_scale,
               pr_gt0 = mean(d > 0), basis = bas, stringsAsFactors = FALSE)
  }
  if (all_contrasts) {
    pairs <- expand.grid(t1 = trts, t2 = trts, stringsAsFactors = FALSE)
    pairs <- pairs[pairs$t1 != pairs$t2, ]
    out <- do.call(rbind, Map(build, pairs$t1, pairs$t2))
  } else {
    others <- setdiff(trts, reference)
    out <- do.call(rbind, lapply(others, build, t2 = reference))
  }
  rownames(out) <- NULL
  attr(out, "sm") <- object$sm
  attr(out, "backtransf") <- backtransf && logsm
  attr(out, "scale") <- output_scale
  attr(out, "target") <- x
  attr(out, "target_mean") <- x
  attr(out, "estimand") <- "average_conditional_link"
  attr(out, "diagnostic_status") <- object$diagnostics$status %||% "unknown"
  class(out) <- c("cpaic_effects", "data.frame")
  out
}

#' @export
relative_effects.cpaic_bridge <- function(object, reference = NULL,
                                          all_contrasts = FALSE,
                                          backtransf = TRUE, level = 0.95,
                                          newdata = NULL, ...) {
  fit <- object$fit
  suffix <- if (object$effect == "random") "random" else "common"
  TE <- fit[[paste0("TE.", suffix)]]
  seTE <- fit[[paste0("seTE.", suffix)]]
  trts <- rownames(TE)
  if (is.null(reference)) reference <- object$reference
  .cpaic_check_ref_level(reference, trts, level)
  z <- stats::qnorm(1 - (1 - level) / 2)
  sm <- object$sm
  logsm <- .is_log_sm(sm)
  output_scale <- if (backtransf && logsm) "natural" else "link"
  C <- object$connectivity$C
  N <- object$connectivity$null_space

  build <- function(t1, t2) {
    ok <- .cpaic_in_rowspace(matrix(C[t1, ] - C[t2, ], nrow = 1L), N)
    if (!ok) {
      return(data.frame(treatment = t1, comparator = t2, estimate = NA_real_,
                        estimate_link = NA_real_, se_link = NA_real_,
                        lower = NA_real_, upper = NA_real_, scale = output_scale,
                        z = NA_real_, p = NA_real_, stringsAsFactors = FALSE))
    }
    est_link <- TE[t1, t2]
    se <- seTE[t1, t2]
    lo <- est_link - z * se
    hi <- est_link + z * se
    zval <- est_link / se
    p <- 2 * stats::pnorm(-abs(zval))
    est <- est_link
    if (backtransf && logsm) {
      est <- exp(est); lo <- exp(lo); hi <- exp(hi)
    }
    data.frame(treatment = t1, comparator = t2, estimate = est,
               estimate_link = est_link, se_link = se,
               lower = lo, upper = hi, scale = output_scale, z = zval, p = p,
               stringsAsFactors = FALSE)
  }

  if (all_contrasts) {
    pairs <- expand.grid(t1 = trts, t2 = trts, stringsAsFactors = FALSE)
    pairs <- pairs[pairs$t1 != pairs$t2, ]
    out <- do.call(rbind, Map(build, pairs$t1, pairs$t2))
  } else {
    others <- setdiff(trts, reference)
    out <- do.call(rbind, lapply(others, build, t2 = reference))
  }
  rownames(out) <- NULL
  attr(out, "sm") <- sm
  attr(out, "backtransf") <- backtransf && logsm
  attr(out, "scale") <- output_scale
  attr(out, "estimand") <- "input_contrast_scale"
  attr(out, "bridge_validity") <- object$bridge_validity %||% NULL
  class(out) <- c("cpaic_effects", "data.frame")
  out
}

#' @export
print.cpaic_effects <- function(x, digits = 3, ...) {
  sm <- attr(x, "sm")
  if (is.null(sm)) sm <- ""
  bt <- isTRUE(attr(x, "backtransf"))
  tgt <- attr(x, "target_mean") %||% attr(x, "target")
  cat("Relative effects (", sm, if (bt) ", back-transformed" else
      ", link scale", ")\n", sep = "")
  if (!is.null(tgt) && length(tgt)) {
    cat("  Average conditional link-scale effect at target means: ",
        paste(names(tgt), signif(tgt, 3), sep = " = ", collapse = ", "),
        "\n", sep = "")
  }
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE)
  if (anyNA(df$estimate)) {
    cat("  NA = not uniquely estimable from this component design",
        " (see estimable_effects()).\n", sep = "")
  }
  if (bt) {
    cat("  `se_link` is on the link (log) scale; the interval is ",
        "back-transformed.\n", sep = "")
  }
  if (!is.null(df$basis) &&
      any(df$basis == "first-order screen", na.rm = TRUE)) {
    cat("  basis \"first-order screen\" = estimable by the row-space criterion",
        " but leaning\n  on aggregate arms or a survival baseline, so it can be",
        " optimistic; check\n  with prior_sensitivity() / estimable_effects_at().\n",
        sep = "")
  }
  invisible(x)
}

#' League table of all pairwise relative effects
#'
#' @param object A `cpaic_bridge` / `cpaic_fit` object.
#' @param backtransf,level See [relative_effects()].
#' @param digits Rounding for the printed cells.
#' @param ... Passed to [relative_effects()] (e.g. `newdata` for [cmlnmr()]).
#' @return A character matrix (treatments x treatments); cell `[i, j]` is
#'   the effect of the row treatment versus the column treatment with its
#'   confidence interval. Non-estimable cells are empty.
#' @export
league_table <- function(object, backtransf = TRUE, level = 0.95,
                         digits = 2, ...) {
  re <- relative_effects(object, all_contrasts = TRUE, backtransf = backtransf,
                         level = level, ...)
  trts <- sort(unique(c(re$treatment, re$comparator)))
  M <- matrix("", length(trts), length(trts), dimnames = list(trts, trts))
  diag(M) <- trts
  for (k in seq_len(nrow(re))) {
    if (is.na(re$estimate[k])) next
    i <- re$treatment[k]; j <- re$comparator[k]
    M[i, j] <- sprintf("%.*f (%.*f, %.*f)", digits, re$estimate[k],
                       digits, re$lower[k], digits, re$upper[k])
  }
  M
}
