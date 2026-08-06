# Target-mean treatment and component hierarchies ------------------------------
#
# Wigle et al. (2026) set out a workflow for building treatment and component
# hierarchies in an (aggregate-data) component NMA:
#
#   Step 1  state the idealized hierarchy question, i.e. the set S to be ranked
#           and the criterion (the ranking metric);
#   Step 2  determine which of the required relative effects are estimable;
#   Step 3  refine S to the estimable subset, or decline to rank;
#   Step 4  compute the ranking metrics and report them alongside the effects.
#
# The reported component contrast is evaluated at effect-modifier means,
#
#     theta_t(x) = C_t' (beta + Gamma x),
#
# so the ranking changes with those means. This is a hierarchy of average
# conditional link-scale effects. It is not a hierarchy of marginal standardized
# ORs, RRs, or HRs.
#
# Step 2 acquires a target-specific analogue as well. The estimable set is a
# function of the target means (see R/estimability.R), so the set that can
# legitimately be ranked can differ between targets. This function performs
# Steps 2 to 4 automatically and reports what it dropped.

.cpaic_rank_engine <- function(object, newdata, what, set,
                               lower_is_better, include_screen_only) {
  stopifnot(inherits(object, "cpaic_mlnmr"))
  C <- object$C.matrix
  ems <- object$effect_modifiers
  x <- .cpaic_target_x(newdata, ems, object$margins)
  Beff <- .cpaic_beta_at(object, x)

  N <- .cpaic_null_space(object$joint_design)
  if (what == "treatment") {
    ref <- object$reference
    elems <- rownames(C)
    Theta <- Beff %*% t(C)
    colnames(Theta) <- rownames(C)
    Draws <- Theta - Theta[, ref]
    Lmat <- C - matrix(C[ref, ], nrow = nrow(C), ncol = ncol(C), byrow = TRUE)
    rownames(Lmat) <- rownames(C)
  } else {
    elems <- object$comps
    Draws <- Beff
    colnames(Draws) <- elems
    Lmat <- diag(ncol(C))
    rownames(Lmat) <- elems
  }

  V <- do.call(rbind, lapply(seq_along(elems), function(i) {
    .cpaic_target_vec(Lmat[i, ], x)
  }))
  ok <- .cpaic_in_rowspace(V, N)
  by_ipd <- .cpaic_in_rowspace(V, .cpaic_null_space(object$joint_design_ipd))
  names(ok) <- names(by_ipd) <- elems

  if (!is.null(set)) {
    bad <- setdiff(set, elems)
    if (length(bad)) {
      stop("`set` names elements not in the network: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    keep0 <- elems %in% set
    elems <- elems[keep0]; ok <- ok[keep0]; by_ipd <- by_ipd[keep0]
    Draws <- Draws[, elems, drop = FALSE]
  }

  dropped <- elems[!ok]
  if (length(dropped)) {
    warning("Dropped from the hierarchy as not estimable at these target ",
            "means: ", paste(dropped, collapse = ", "),
            ". Ranking them would rank the prior. See estimable_effects_at().",
            call. = FALSE)
  }
  dropped_screen <- character(0)
  keep <- ok
  if (!include_screen_only) {
    dropped_screen <- elems[ok & !by_ipd]
    if (length(dropped_screen)) {
      warning("Dropped from the hierarchy as identified only by aggregate arms ",
              "(a first-order screen that can be optimistic): ",
              paste(dropped_screen, collapse = ", "),
              ". Set include_screen_only = TRUE to rank them as an explicitly ",
              "exploratory hierarchy.", call. = FALSE)
    }
    keep <- ok & by_ipd
  }
  elems <- elems[keep]
  if (length(elems) < 2L) {
    stop("Fewer than two elements are estimable",
         if (!include_screen_only)
           " (excluding elements identified only by aggregate arms)" else "",
         " at these target means, so no hierarchy can be formed. See ",
         "estimable_effects_at()",
         if (!include_screen_only) " or set include_screen_only = TRUE" else "",
         ".", call. = FALSE)
  }
  Draws <- Draws[, elems, drop = FALSE]

  sgn <- if (lower_is_better) 1 else -1
  R <- t(apply(sgn * Draws, 1L, rank, ties.method = "average"))
  n <- length(elems)
  P <- matrix(0, nrow = n, ncol = n,
              dimnames = list(elems, as.character(seq_len(n))))
  for (i in seq_len(nrow(Draws))) {
    ord <- order(sgn * Draws[i, ])
    groups <- split(seq_len(n),
                    cumsum(c(TRUE, diff((sgn * Draws[i, ])[ord]) != 0)))
    for (positions in groups) {
      P[ord[positions], positions] <-
        P[ord[positions], positions] + 1 / length(positions)
    }
  }
  P <- P / nrow(Draws)

  list(
    elements = elems,
    draws = Draws,
    ranks = R,
    probabilities = P,
    target = x,
    dropped = dropped,
    dropped_screen = dropped_screen
  )
}

#' Treatment and component hierarchies at target covariate means
#'
#' Ranks treatments or components using average conditional link-scale effects
#' evaluated at supplied effect-modifier means. Because the component contrast
#' is linear in those means, the hierarchy can change with them. This is not a
#' hierarchy of marginal effects standardized over a covariate distribution.
#'
#' Elements whose relative effect is not estimable at those target means are
#' **dropped from the ranking set** rather than ranked from a prior-driven
#' posterior, and are reported in the `dropped` attribute. This is Step 3 of the
#' Wigle et al. workflow, and it matters more here than in the aggregate-data
#' case, because the estimable set depends on the target means (see
#' [estimable_effects_at()]).
#'
#' Ranking metrics depend on the set being ranked, so they are not comparable
#' across different sets. Report them alongside the relative effects, never
#' instead of them.
#'
#' @param object A [cmlnmr()] fit.
#' @param newdata A one-row data frame giving target effect-modifier means.
#'   Required when the model has effect modifiers.
#' @param what `"treatment"` (default) or `"component"`. Ranking components by
#'   their incremental effect is only meaningful in an additive model.
#' @param set Optional character vector restricting the elements to rank (the
#'   set `S` of Wigle et al.). Defaults to all treatments (or all components).
#' @param lower_is_better If `TRUE`, a smaller effect is preferred (e.g.
#'   mortality). Default `FALSE` (a larger effect is preferred).
#' @param include_screen_only If `FALSE` (default), elements whose relative
#'   effect is identified only by aggregate arms (a first-order screen that can
#'   be optimistic under a nonlinear link) are excluded from the hierarchy and
#'   reported in the `dropped_screen` attribute. Set `TRUE` to rank them as an
#'   explicitly exploratory hierarchy.
#'
#'   The test here is whether the individual patient data identify the element,
#'   which is not identical to `basis == "exact"` in [estimable_effects_at()].
#'   That column additionally excludes **survival** from `"exact"`, because a
#'   flexible baseline hazard adds support-dependent nuisance parameters the
#'   covariate-support argument does not see. Survival elements identified by IPD
#'   are therefore still ranked by default; dropping every survival element from
#'   every survival hierarchy would leave nothing to rank. Read a survival
#'   hierarchy alongside [estimable_effects_at()] rather than on its own.
#' @param estimand The only implemented value is
#'   `"average_conditional_link"`. Marginal standardized rankings are not yet
#'   implemented and are rejected explicitly.
#' @param ... Unused.
#'
#' @return A data frame, ordered from most to least preferred, with columns
#'   `element`, `estimate` (posterior mean of the relative effect versus the
#'   reference, on the link scale), `p_best`, `median_rank`, `mean_rank` and
#'   `sucra`. The `dropped` attribute lists elements excluded as not estimable
#'   at these target means.
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @seealso [estimable_effects_at()], [relative_effects()]
#' @examplesIf FALSE
#' # Which component ranks best when the target mean of x1 is 0.5?
#' cpaic_ranks(fit, newdata = data.frame(x1 = 0.5), what = "component")
#' @export
cpaic_ranks <- function(object, newdata = NULL,
                        what = c("treatment", "component"), set = NULL,
                        lower_is_better = FALSE, include_screen_only = FALSE,
                        estimand = "average_conditional_link", ...) {
  .cpaic_match_estimand(estimand, "average_conditional_link", "cpaic_ranks()")
  what <- match.arg(what)
  ranked <- .cpaic_rank_engine(
    object = object, newdata = newdata, what = what, set = set,
    lower_is_better = lower_is_better,
    include_screen_only = include_screen_only
  )
  n <- length(ranked$elements)
  mean_rank <- as.numeric(ranked$probabilities %*% seq_len(n))
  out <- data.frame(
    element = ranked$elements,
    estimate = colMeans(ranked$draws),
    p_best = ranked$probabilities[, 1L],
    median_rank = apply(ranked$ranks, 2L, stats::median),
    mean_rank = mean_rank,
    sucra = (n - mean_rank) / (n - 1),
    row.names = NULL, stringsAsFactors = FALSE
  )
  out <- out[order(out$mean_rank), ]
  rownames(out) <- NULL
  attr(out, "dropped") <- ranked$dropped
  attr(out, "dropped_screen") <- ranked$dropped_screen
  attr(out, "target") <- ranked$target
  attr(out, "target_mean") <- ranked$target
  attr(out, "estimand") <- "average_conditional_link"
  attr(out, "diagnostic_status") <- object$diagnostics$status %||% "unknown"
  attr(out, "what") <- what
  class(out) <- c("cpaic_ranks", "data.frame")
  out
}

#' @export
print.cpaic_ranks <- function(x, digits = 3, ...) {
  tgt <- attr(x, "target_mean") %||% attr(x, "target")
  cat("Average conditional link-scale ", attr(x, "what"),
      " hierarchy\n", sep = "")
  if (length(tgt)) {
    cat("  Target effect-modifier means: ",
        paste(names(tgt), signif(tgt, 3), sep = " = ", collapse = ", "),
        "\n", sep = "")
  }
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE)
  dr <- attr(x, "dropped")
  if (length(dr)) {
    cat("  Not estimable at these means, so not ranked: ",
        paste(dr, collapse = ", "), "\n", sep = "")
  }
  ds <- attr(x, "dropped_screen")
  if (length(ds)) {
    cat("  Identified only by aggregate arms (first-order screen), so not ",
        "ranked by default: ", paste(ds, collapse = ", "),
        "\n  (set include_screen_only = TRUE to rank them as exploratory).\n",
        sep = "")
  }
  cat("  Ranking metrics depend on the set ranked; report them with the",
      " effects, not instead.\n", sep = "")
  invisible(x)
}

#' How a hierarchy changes across target effect-modifier means
#'
#' Recomputes [cpaic_ranks()] over a grid of target means. This exposes how the
#' hierarchy of average conditional link-scale effects changes with the chosen
#' mean. It does not standardize effects over a sequence of target distributions.
#'
#' @param object A [cmlnmr()] fit.
#' @param em Name of the effect modifier to vary.
#' @param values Numeric vector of target values for `em`.
#' @param at Optional named vector fixing the other effect modifiers. Defaults
#'   to 0 for each.
#' @param what,lower_is_better,include_screen_only See [cpaic_ranks()].
#' @param ... Unused.
#' @return A data frame with one row per (element, target value), giving `sucra`,
#'   `mean_rank`, `p_best`, and `estimate`. Failed target values are retained as
#'   one `status = "failed"` row with `NA` metrics and an explanatory `error`.
#' @seealso [cpaic_ranks()]
#' @examplesIf FALSE
#' rank_curve(fit, em = "x1", values = seq(-1, 1, by = 0.25), what = "component")
#' @export
rank_curve <- function(object, em, values, at = NULL,
                       what = c("treatment", "component"),
                       lower_is_better = FALSE, include_screen_only = FALSE,
                       ...) {
  stopifnot(inherits(object, "cpaic_mlnmr"))
  what <- match.arg(what)
  ems <- object$effect_modifiers
  if (!em %in% ems) {
    stop("`em` must be one of the effect modifiers: ",
         paste(ems, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(values) || !length(values) || any(!is.finite(values))) {
    stop("`values` must be a non-empty numeric vector of finite target values.",
         call. = FALSE)
  }
  if (!is.null(at)) {
    if (!is.numeric(at) || is.null(names(at)) || any(!nzchar(names(at))) ||
        anyDuplicated(names(at)) || any(!is.finite(at))) {
      stop("`at` must be a named numeric vector of finite target values.",
           call. = FALSE)
    }
    bad <- setdiff(names(at), ems)
    if (length(bad)) {
      stop("`at` names effect modifiers not in the model: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    if (em %in% names(at)) {
      stop("`at` must not include the effect modifier varied by `values`.",
           call. = FALSE)
    }
  }
  base <- stats::setNames(rep(0, length(ems)), ems)
  if (!is.null(at)) base[names(at)] <- at

  res <- lapply(seq_along(values), function(i) {
    v <- values[i]
    nd <- as.data.frame(as.list(replace(base, em, v)))
    r <- tryCatch(
      suppressWarnings(cpaic_ranks(object, newdata = nd, what = what,
                                   lower_is_better = lower_is_better,
                                   include_screen_only = include_screen_only)),
      error = identity)
    if (inherits(r, "error")) {
      return(data.frame(element = NA_character_, value = v,
                        estimate = NA_real_, p_best = NA_real_,
                        mean_rank = NA_real_, sucra = NA_real_,
                        status = "failed", error = conditionMessage(r),
                        target_index = i, row.names = NULL,
                        stringsAsFactors = FALSE))
    }
    data.frame(element = r$element, value = v, estimate = r$estimate,
               p_best = r$p_best, mean_rank = r$mean_rank, sucra = r$sucra,
               status = "ok", error = NA_character_, target_index = i,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  names(out)[names(out) == "value"] <- em
  attr(out, "em") <- em
  attr(out, "what") <- what
  attr(out, "estimand") <- "average_conditional_link"
  out
}
