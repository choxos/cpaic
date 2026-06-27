# Post-fit relative effects and league tables ---------------------------------

#' Is the summary measure on a (natural) log scale?
#' @noRd
.is_log_sm <- function(sm) {
  toupper(sm) %in% c("OR", "RR", "HR", "IRR", "RRR", "ROR", "HRR")
}

#' Relative treatment effects from a cpaic fit
#'
#' Tidies the (random- or common-effects) relative effects from the
#' component-bridged model: every treatment versus a chosen reference, or
#' all pairwise comparisons. Effects are reported on the natural scale of
#' the summary measure (e.g. odds ratios) unless `backtransf = FALSE`.
#'
#' @param object A `cpaic_bridge`, `cpaic_maic`, or `cpaic_stc` object.
#' @param reference Reference treatment. Defaults to the network reference.
#' @param all_contrasts If `TRUE`, return all pairwise comparisons instead
#'   of versus the reference.
#' @param backtransf If `TRUE` (default) back-transform log-scale measures
#'   (OR/RR/HR/...) by exponentiating.
#' @param level Confidence level for the intervals. Default `0.95`.
#' @param ... Unused.
#'
#' @return A data frame with columns `treatment`, `comparator`, `estimate`,
#'   `se` (link scale), `lower`, `upper`, `z`, `p`.
#' @export
relative_effects <- function(object, reference = NULL, all_contrasts = FALSE,
                             backtransf = TRUE, level = 0.95, ...) {
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
                                         ...) {
  C <- object$C.matrix
  beta <- as.matrix(object$fit$draws("beta", format = "draws_matrix"))
  Theta <- beta %*% t(C)              # draws x treatments (theta = C beta)
  colnames(Theta) <- rownames(C)
  trts <- rownames(C)
  if (is.null(reference)) reference <- object$reference
  a <- (1 - level) / 2
  logsm <- .is_log_sm(object$sm)
  bt <- function(v) if (backtransf && logsm) exp(v) else v

  build <- function(t1, t2) {
    d <- Theta[, t1] - Theta[, t2]
    q <- stats::quantile(d, c(a, 1 - a), names = FALSE)
    data.frame(treatment = t1, comparator = t2, estimate = bt(mean(d)),
               se = stats::sd(d), lower = bt(q[1]), upper = bt(q[2]),
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
  class(out) <- c("cpaic_effects", "data.frame")
  out
}

#' @export
relative_effects.cpaic_bridge <- function(object, reference = NULL,
                                          all_contrasts = FALSE,
                                          backtransf = TRUE, level = 0.95,
                                          ...) {
  fit <- object$fit
  suffix <- if (object$effect == "random") "random" else "common"
  TE <- fit[[paste0("TE.", suffix)]]
  seTE <- fit[[paste0("seTE.", suffix)]]
  trts <- rownames(TE)
  if (is.null(reference)) reference <- object$reference
  z <- stats::qnorm(1 - (1 - level) / 2)
  sm <- object$sm
  logsm <- .is_log_sm(sm)

  build <- function(t1, t2) {
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
  bt <- attr(x, "backtransf")
  cat("Relative effects (", sm, if (bt) ", back-transformed" else
      ", link scale", ")\n", sep = "")
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE)
  invisible(x)
}

#' League table of all pairwise relative effects
#'
#' @param object A `cpaic_bridge` / `cpaic_fit` object.
#' @param backtransf,level See [relative_effects()].
#' @param digits Rounding for the printed cells.
#' @return A character matrix (treatments x treatments); cell `[i, j]` is
#'   the effect of the row treatment versus the column treatment with its
#'   confidence interval.
#' @export
league_table <- function(object, backtransf = TRUE, level = 0.95,
                         digits = 2) {
  re <- relative_effects(object, all_contrasts = TRUE, backtransf = backtransf,
                         level = level)
  trts <- sort(unique(c(re$treatment, re$comparator)))
  M <- matrix("", length(trts), length(trts), dimnames = list(trts, trts))
  diag(M) <- trts
  for (k in seq_len(nrow(re))) {
    i <- re$treatment[k]; j <- re$comparator[k]
    M[i, j] <- sprintf("%.*f (%.*f, %.*f)", digits, re$estimate[k],
                       digits, re$lower[k], digits, re$upper[k])
  }
  M
}
