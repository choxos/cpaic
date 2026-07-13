# Estimability of population-adjusted component contrasts ---------------------
#
# Wigle et al. (2026) show that in an aggregate-data component NMA the uniquely
# estimable relative effects are exactly the row space of the design matrix
# X = B C. Population adjustment changes the question, because the estimand
#
#     theta_t(x) - theta_u(x) = m' (beta + Gamma x),      m = C_t - C_u,
#
# involves the component x effect-modifier interactions Gamma as well as the
# component main effects beta. Estimability of m'beta is then necessary but not
# sufficient for the population-adjusted contrast at any x != 0.
#
# Stack the parameters as phi = (beta, vec(Gamma)) in R^{K(1+Q)}. The estimand is
# the linear functional v' phi with
#
#     v(m, x) = (1, x') %x% m                                            (*)
#
# so estimability is again a row-space question, now in the augmented parameter
# space. The first-order information design D is built from:
#
#   * each two-arm IPD study s (contrast m_s): the row (e_0 %x% m_s), giving
#     m_s'beta, plus one row ((0, a') %x% m_s) for every effect-modifier
#     direction a that actually VARIES within that study, giving m_s'Gamma a.
#     Within an IPD study you can regress on arm and arm-by-covariate, so the
#     arm main effect and each arm-by-covariate slope are identified directly.
#
#   * each aggregate contrast s (contrast m_s, covariate mean xbar_s): the single
#     row ((1, xbar_s') %x% m_s), giving m_s'(beta + Gamma xbar_s), which is the
#     contrast in that study's OWN population.
#
# A contrast is then identified by the first-order information iff
# v(m, x) is in rowspace(D).
#
# Properties (all checked in tests/testthat/test-estimability.R):
#   Q = 0 collapses to the Wigle et al. (2026) criterion exactly.
#   For an identity link the criterion is exact (the aggregate contrast really is
#     m_s'(beta + Gamma xbar_s)).
#   For a nonlinear link (logit, log) the criterion is a DESIGN-BASED SCREEN, not
#     a complete identification theorem, and it is NOT correct to call it
#     "conservative". The aggregate arm likelihood is a nonlinear integral, so the
#     exact local information is the Jacobian of the mean map, in which an
#     aggregate study pins the contrast down at a VARIANCE-WEIGHTED mean rather
#     than at the raw mean xbar_s. The screen has the right rank structure (a
#     two-arm aggregate study supplies exactly one functional per contrast
#     direction, which cannot separate the 1 + Q unknowns), so it finds the
#     under-determined contrasts correctly; but the anchor point is shifted. Do
#     not read a failure as "weakly identified through curvature": the
#     nonlinearity generally identifies other, NONLINEAR functionals of phi (such
#     as the marginal natural-scale effect), not the conditional link-scale
#     contrast m'(beta + Gamma x).
#   An observed aggregate contrast is always estimable in its own population
#     (take x = xbar_s in (*)).
#   A two-arm IPD study in which every effect modifier varies identifies its own
#     contrast in EVERY target population, which is what anchored STC/MAIC do.
#   Estimability can depend on the target population: the estimable set can
#     shrink as x moves away from the covariate origin.

#' Augmented target-contrast vector v = (1, x) %x% m
#' @noRd
.cpaic_target_vec <- function(m, x) {
  as.numeric(kronecker(matrix(c(1, x), nrow = 1L), matrix(m, nrow = 1L)))
}

#' First-order information design for (beta, vec(Gamma))
#'
#' @param C treatment-by-component matrix.
#' @param ipd data frame of IPD (already validated).
#' @param agd data frame of aggregate arms.
#' @param effect_modifiers character vector of EM names.
#' @param study,trt column names.
#' @param tol tolerance for detecting within-study covariate variation.
#' @return A matrix with `ncol(C) * (1 + Q)` columns.
#' @noRd
.cpaic_joint_design <- function(C, ipd, agd, effect_modifiers,
                                study = ".study", trt = ".trt", tol = 1e-8) {
  K <- ncol(C)
  Q <- length(effect_modifiers)
  rows <- list()

  contrasts_of <- function(ts) {
    # All arm contrasts of a study, relative to its first arm.
    Cs <- C[match(ts, rownames(C)), , drop = FALSE]
    sweep(Cs[-1L, , drop = FALSE], 2, Cs[1L, ], "-")
  }

  # --- IPD studies: arm main effect + arm-by-covariate slopes ----------------
  if (!is.null(ipd) && nrow(ipd)) {
    for (ss in unique(as.character(ipd[[study]]))) {
      sub <- ipd[as.character(ipd[[study]]) == ss, , drop = FALSE]
      ts <- unique(as.character(sub[[trt]]))
      if (length(ts) < 2L) next
      M <- contrasts_of(ts)                       # (n_arms - 1) x K

      # The LEVEL row must be anchored at this study's OWN covariate mean, not at
      # the global covariate origin. If a modifier is held constant within the
      # trial (a single-sex trial, say), the trial identifies the arm effect only
      # AT that value: with x == 5 one learns m'(beta + 5 Gamma) and can never
      # separate m'beta from m'Gamma. Anchoring at the origin would call x = 0
      # estimable (it is not) and x = 5 not estimable (it is), i.e. exactly
      # backwards. When every modifier varies, the level and slope rows together
      # span the full block I_{1+Q} %x% m, so nothing is lost.
      xbar <- rep(0, Q)
      A <- matrix(numeric(0), nrow = 0L, ncol = Q)
      if (Q) {
        Xs <- as.matrix(sub[, effect_modifiers, drop = FALSE])
        xbar <- colMeans(Xs)
        Xc <- sweep(Xs, 2, xbar, "-")
        sv <- svd(Xc, nu = 0L, nv = Q)
        keep <- sv$d > tol * max(1, max(sv$d))
        if (any(keep)) A <- t(sv$v[, keep, drop = FALSE])   # r x Q basis
      }

      for (r in seq_len(nrow(M))) {
        m <- M[r, ]
        # level, anchored at this study's own covariate mean
        rows[[length(rows) + 1L]] <- .cpaic_target_vec(m, xbar)
        # slopes: (0, a) %x% m for each direction that actually varies
        for (j in seq_len(nrow(A))) {
          a <- A[j, ]
          rows[[length(rows) + 1L]] <-
            as.numeric(kronecker(matrix(c(0, a), nrow = 1L),
                                 matrix(m, nrow = 1L)))
        }
      }
    }
  }

  # --- Aggregate studies: one row per contrast, at that study's own mean -----
  if (!is.null(agd) && nrow(agd)) {
    for (ss in unique(as.character(agd[[study]]))) {
      sub <- agd[as.character(agd[[study]]) == ss, , drop = FALSE]
      ts <- as.character(sub[[trt]])
      if (length(ts) < 2L) next
      M <- contrasts_of(ts)
      xbar <- if (Q) {
        vapply(effect_modifiers,
               function(v) mean(sub[[paste0(v, "_mean")]]), numeric(1))
      } else numeric(0)
      for (r in seq_len(nrow(M))) {
        rows[[length(rows) + 1L]] <- .cpaic_target_vec(M[r, ], xbar)
      }
    }
  }

  if (!length(rows)) return(matrix(0, nrow = 1L, ncol = K * (1L + Q)))
  do.call(rbind, rows)
}

#' Which population-adjusted contrasts are estimable at a target population?
#'
#' Extends the row-space criterion of Wigle et al. (2026) from the component
#' main effects to the population-adjusted estimand
#' `theta_t(x) = C_t' (beta + Gamma x)`. A relative effect is identified by the
#' first-order information if and only if its augmented contrast vector
#' `(1, x) %x% (C_t - C_u)` lies in the row space of the information design
#' (see the file header for how that design is built from the IPD and aggregate
#' evidence).
#'
#' Because the criterion depends on `x`, **the estimable set can depend on the
#' target population**: a contrast estimable at the covariate origin need not be
#' estimable in a target population where the component by effect-modifier
#' interactions are not identified.
#'
#' @param object A `cpaic_mlnmr` fit.
#' @param newdata A one-row data frame giving the target population's
#'   effect-modifier values. Defaults to the covariate origin.
#' @param reference Reference treatment. Defaults to the fit's reference.
#' @param ... Unused.
#' @return A data frame with `treatment`, `comparator`, `estimable`, and
#'   `identified_by` (`"IPD"`, `"aggregate"`, or `"none"`).
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @export
estimable_effects_at <- function(object, newdata = NULL, reference = NULL,
                                 ...) {
  stopifnot(inherits(object, "cpaic_mlnmr"))
  C <- object$C.matrix
  ems <- object$effect_modifiers
  Q <- length(ems)
  if (is.null(reference)) reference <- object$reference
  x <- if (is.null(newdata)) rep(0, Q) else
    .cpaic_target_x(newdata, ems)

  D <- object$joint_design
  N <- .cpaic_null_space(D)
  trts <- rownames(C)
  others <- setdiff(trts, reference)

  V <- do.call(rbind, lapply(others, function(t1) {
    .cpaic_target_vec(C[t1, ] - C[reference, ], x)
  }))
  ok <- .cpaic_in_rowspace(V, N)

  # Was the interaction information supplied by IPD or only by aggregate arms?
  D_ipd <- object$joint_design_ipd
  N_ipd <- .cpaic_null_space(D_ipd)
  by_ipd <- .cpaic_in_rowspace(V, N_ipd)

  data.frame(
    treatment = others,
    comparator = reference,
    estimable = ok,
    identified_by = ifelse(!ok, "none", ifelse(by_ipd, "IPD", "aggregate")),
    row.names = NULL, stringsAsFactors = FALSE
  )
}
