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

#' Resolve a target population into a vector of effect-modifier values
#'
#' ML-NMR relative effects are population-specific whenever effect modifiers
#' are present: there is no population-free relative effect. `newdata` must
#' therefore name the target population's effect-modifier values.
#' @noRd
.cpaic_target_x <- function(newdata, effect_modifiers) {
  Q <- length(effect_modifiers)
  if (!Q) return(numeric(0))
  if (is.null(newdata)) {
    stop("Specify `newdata`: a one-row data frame giving the target ",
         "population's effect-modifier values (",
         paste(effect_modifiers, collapse = ", "), "). Relative effects from ",
         "a model with effect modifiers are population-specific, so there is ",
         "no population-free answer. Use, for example, ",
         "newdata = data.frame(", effect_modifiers[1], " = 0) to obtain the ",
         "effect at the covariate origin.", call. = FALSE)
  }
  newdata <- as.data.frame(newdata)
  if (nrow(newdata) != 1L) {
    stop("`newdata` must have exactly one row (the target population).",
         call. = FALSE)
  }
  miss <- setdiff(effect_modifiers, names(newdata))
  if (length(miss)) {
    stop("`newdata` is missing effect modifier(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  x <- vapply(effect_modifiers, function(v) as.numeric(newdata[[v]][1]),
              numeric(1))
  if (any(!is.finite(x))) {
    stop("`newdata` effect-modifier values must be finite.", call. = FALSE)
  }
  x
}

#' Posterior draws of the component effects in a target population
#'
#' In the component-additive ML-NMR the treatment effect of treatment `t` in a
#' population with effect-modifier values `x` is
#' `theta_t(x) = C_t' (beta + gamma x)`, so `beta + gamma x` are the component
#' effects *in that population*. Reporting `beta` alone (the effect at the
#' covariate origin) is not a population-adjusted quantity.
#' @noRd
.cpaic_beta_at <- function(object, x) {
  B <- as.matrix(object$fit$draws("beta", format = "draws_matrix"))
  nC <- ncol(object$C.matrix)
  if (!length(x)) return(B)
  G <- as.matrix(object$fit$draws("gamma", format = "draws_matrix"))
  out <- B
  for (q in seq_along(x)) {
    cols <- paste0("gamma[", seq_len(nC), ",", q, "]")
    miss <- setdiff(cols, colnames(G))
    if (length(miss)) {
      stop("The fit does not contain the component x effect-modifier draws ",
           "needed to evaluate a target population.", call. = FALSE)
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
#' interactions, so relative effects are **population-specific**:
#' `theta_t(x) = C_t' (beta + gamma x)`. You must name the target population
#' through `newdata`; there is no population-free relative effect.
#'
#' @param object A fitted cpaic object (`cpaic_bridge`, `cpaic_maic`,
#'   `cpaic_stc`, or `cpaic_mlnmr`).
#' @param reference Reference treatment. Defaults to the network reference.
#' @param all_contrasts If `TRUE`, return all pairwise comparisons instead
#'   of versus the reference.
#' @param backtransf If `TRUE` (default) back-transform log-scale measures
#'   (OR/RR/HR/...) by exponentiating.
#' @param level Confidence level for the intervals. Default `0.95`.
#' @param newdata For [cmlnmr()] fits: a one-row data frame giving the target
#'   population's effect-modifier values. Required when the model has effect
#'   modifiers.
#' @param ... Unused.
#'
#' @return A data frame with columns `treatment`, `comparator`, `estimate`,
#'   `se` (link scale), `lower`, `upper`, and `z`/`p` for frequentist fits.
#'   For [cmlnmr()] (Bayesian) fits the intervals are credible intervals and
#'   the final column is `pr_gt0`, the posterior probability that the
#'   effect (on the link scale) exceeds zero, instead of `z`/`p`.
#' @export
relative_effects <- function(object, reference = NULL, all_contrasts = FALSE,
                             backtransf = TRUE, level = 0.95, newdata = NULL,
                             ...) {
  UseMethod("relative_effects")
}

#' @export
relative_effects.cpaic_fit <- function(object, ...) {
  relative_effects(object$bridge, ...)
}

#' @export
relative_effects.cpaic_mlnmr <- function(object, reference = NULL,
                                         all_contrasts = FALSE,
                                         backtransf = TRUE, level = 0.95,
                                         newdata = NULL, ...) {
  C <- object$C.matrix
  trts <- rownames(C)
  if (is.null(reference)) reference <- object$reference
  .cpaic_check_ref_level(reference, trts, level)

  x <- .cpaic_target_x(newdata, object$effect_modifiers)
  Beff <- .cpaic_beta_at(object, x)      # draws x components, in population x
  Theta <- Beff %*% t(C)                 # draws x treatments
  colnames(Theta) <- trts

  a <- (1 - level) / 2
  logsm <- .is_log_sm(object$sm)
  # Estimability of a POPULATION-ADJUSTED contrast depends on the target
  # population: it needs the relevant rows of Gamma, not just of beta. Test the
  # augmented contrast (1, x) %x% m against the joint information design.
  Njoint <- .cpaic_null_space(object$joint_design)

  build <- function(t1, t2) {
    v <- .cpaic_target_vec(C[t1, ] - C[t2, ], x)
    ok <- .cpaic_in_rowspace(matrix(v, nrow = 1L), Njoint)
    if (!ok) {
      return(data.frame(treatment = t1, comparator = t2, estimate = NA_real_,
                        se = NA_real_, lower = NA_real_, upper = NA_real_,
                        pr_gt0 = NA_real_, stringsAsFactors = FALSE))
    }
    d <- Theta[, t1] - Theta[, t2]                 # link scale
    e <- if (backtransf && logsm) exp(d) else d     # reporting scale
    q <- stats::quantile(e, c(a, 1 - a), names = FALSE)
    data.frame(treatment = t1, comparator = t2, estimate = mean(e),
               se = stats::sd(d), lower = q[1], upper = q[2],
               pr_gt0 = mean(d > 0), stringsAsFactors = FALSE)
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
  attr(out, "target") <- x
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
  C <- object$connectivity$C
  N <- object$connectivity$null_space

  build <- function(t1, t2) {
    ok <- .cpaic_in_rowspace(matrix(C[t1, ] - C[t2, ], nrow = 1L), N)
    if (!ok) {
      return(data.frame(treatment = t1, comparator = t2, estimate = NA_real_,
                        se = NA_real_, lower = NA_real_, upper = NA_real_,
                        z = NA_real_, p = NA_real_, stringsAsFactors = FALSE))
    }
    est <- TE[t1, t2]
    se <- seTE[t1, t2]
    lo <- est - z * se
    hi <- est + z * se
    zval <- est / se
    p <- 2 * stats::pnorm(-abs(zval))
    if (backtransf && logsm) {
      est <- exp(est); lo <- exp(lo); hi <- exp(hi)
    }
    data.frame(treatment = t1, comparator = t2, estimate = est, se = se,
               lower = lo, upper = hi, z = zval, p = p,
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
  class(out) <- c("cpaic_effects", "data.frame")
  out
}

#' @export
print.cpaic_effects <- function(x, digits = 3, ...) {
  sm <- attr(x, "sm")
  if (is.null(sm)) sm <- ""
  bt <- isTRUE(attr(x, "backtransf"))
  tgt <- attr(x, "target")
  cat("Relative effects (", sm, if (bt) ", back-transformed" else
      ", link scale", ")\n", sep = "")
  if (!is.null(tgt) && length(tgt)) {
    cat("  Target population: ",
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
