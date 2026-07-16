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
#     under-determined contrasts correctly; but the anchor point is shifted.
#
#     Concretely, an aggregate study need not identify its own-population contrast
#     under a nonlinear link. With a log link, one aggregate study, and a symmetric
#     covariate P(x = -1) = P(x = +1) = 1/2 (so xbar = 0), the arm means are
#     exp(mu) and exp(mu + beta) cosh(gamma), so the data identify only
#     beta + log cosh(gamma), NOT beta = m'(beta + Gamma xbar). Do not read a
#     failure of the screen as "weakly identified through curvature" either: the
#     nonlinearity generally identifies other, NONLINEAR functionals of phi (such
#     as the marginal natural-scale effect), not the conditional link-scale
#     contrast m'(beta + Gamma x).
#
#   ESTIMAND. m'(beta + Gamma x) is the CONDITIONAL contrast on the linear
#     predictor scale, evaluated at covariate value x. It is not the marginal
#     effect obtained by integrating natural-scale outcomes over the target
#     population; on a non-collapsible scale those differ. ML-NMR distinguishes
#     them (Phillippo et al. 2020) and so must any statement built on this
#     criterion, including the hierarchies in R/ranks.R.
#
#   ECOLOGICAL vs WITHIN-STUDY INTERACTION. Writing an effect as
#     alpha + gamma_W (x - xbar_s) + gamma_B xbar_s, aggregate contrasts depend
#     only on alpha + gamma_B xbar_s: they carry NO information about the
#     within-study interaction gamma_W, and the model silently imposes
#     gamma_W = gamma_B. Randomization identifies each study's effect but does not
#     randomize covariate means ACROSS studies, so a between-study gradient is
#     confounded in a way a within-study slope is not (Berlin et al. 2002;
#     Freeman et al. 2018). `identified_by` therefore separates the two.
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

#' Orthonormal basis (columns) of the row space of a matrix
#'
#' For an n x p matrix this returns a basis of the space spanned by its ROWS,
#' living in R^p. (The column space would live in R^n, which is the wrong space
#' for a covariate-direction argument.)
#' @noRd
.cpaic_row_space <- function(X, tol = 1e-8) {
  X <- as.matrix(X)
  if (!nrow(X) || !ncol(X)) return(matrix(0, ncol(X), 0L))
  s <- svd(X, nu = 0L, nv = ncol(X))
  # `s$d` has length min(nrow, ncol), while `s$v` has ncol(X) columns. When
  # nrow < ncol the singular values run out before the columns of `s$v`, so a
  # logical `keep` of length min(nrow, ncol) would RECYCLE across the extra
  # columns and return null-space directions as if they spanned the row space:
  # a thin arm (fewer patients than 1 + Q) would then look fully supported and a
  # prior-only contrast would be called identified. Indexing by position keeps
  # only genuine singular directions and never reaches the padding columns.
  keep <- which(s$d > tol * max(1, max(s$d)))
  if (!length(keep)) return(matrix(0, ncol(X), 0L))
  s$v[, keep, drop = FALSE]
}

#' Orthonormal basis of the intersection of two subspaces
#'
#' `B1`, `B2` have orthonormal columns spanning subspaces of the same R^p.
#' A vector lies in both iff it is `B1 c1 = B2 c2`, i.e. `(c1, c2)` is in the
#' null space of `[B1, -B2]`.
#' @noRd
.cpaic_subspace_intersect <- function(B1, B2, tol = 1e-8) {
  p <- nrow(B1)
  if (!ncol(B1) || !ncol(B2)) return(matrix(0, p, 0L))
  N <- .cpaic_null_space(cbind(B1, -B2), tol = tol)
  if (!ncol(N)) return(matrix(0, p, 0L))
  V <- B1 %*% N[seq_len(ncol(B1)), , drop = FALSE]
  s <- svd(V, nv = 0L)
  keep <- s$d > tol * max(1, max(s$d))
  if (!any(keep)) return(matrix(0, p, 0L))
  s$u[, keep, drop = FALSE]
}

#' Which functionals w'(m'beta, m'Gamma_1, ..., m'Gamma_Q) does an IPD arm
#' contrast identify?
#'
#' The linear predictor carries an arm-COMMON nuisance surface (the study
#' intercept and the prognostic effects), which cancels from an arm contrast
#' evaluated at the same covariate value. What survives is that the contrast is
#' pinned down only on the augmented covariate directions that BOTH arms
#' actually support. Pooling the two arms' covariates would overstate this:
#' if the control arm has x in {0, 1} but the treated arm has x == 5 (perfect
#' confounding of arm with covariate), the pooled sample looks like it varies,
#' yet the only identified functional is m'(beta + 5 Gamma).
#'
#' The identified set is therefore `row(U_ref)` intersect `row(U_arm)`, where
#' `U_a` stacks the augmented rows `(1, x_i)` of arm `a`. In a randomized trial
#' the two arms share a covariate distribution and this equals the pooled row
#' space, so nothing is lost in the usual case.
#' @noRd
.cpaic_arm_directions <- function(Xref, Xarm, tol = 1e-8) {
  aug <- function(X) cbind(1, as.matrix(X))
  B1 <- .cpaic_row_space(aug(Xref), tol = tol)
  B2 <- .cpaic_row_space(aug(Xarm), tol = tol)
  .cpaic_subspace_intersect(B1, B2, tol = tol)
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

  # --- IPD studies ------------------------------------------------------------
  # An arm contrast cancels the arm-common nuisance surface (the study intercept
  # and the prognostic effects) only when evaluated at the same covariate value.
  # So an arm contrast identifies w'(m'beta, m'Gamma_1, ..., m'Gamma_Q) exactly
  # for the augmented covariate directions w that BOTH arms support:
  #
  #     w in row(U_ref) intersect row(U_arm),   U_a = rows (1, x_i) of arm a.
  #
  # Two things follow, and both matter.
  #
  #   * The anchor is the study's own covariate support, never the global origin.
  #     A modifier held constant within a trial (a single-sex trial, say) gives
  #     w = (1, 5): the trial identifies m'(beta + 5 Gamma) and can never
  #     separate m'beta from m'Gamma. Anchoring at the origin would call x = 0
  #     estimable (it is not) and deny x = 5 (which is), i.e. exactly backwards.
  #
  #   * Pooling the arms would overstate identification. If the control arm has
  #     x in {0, 1} but the treated arm has x == 5, the pooled sample appears to
  #     vary, yet arm is perfectly confounded with the covariate and only
  #     m'(beta + 5 Gamma) is identified. In a randomized trial the arms share a
  #     covariate distribution, the intersection equals the pooled row space, and
  #     nothing is lost.
  #
  # For a multi-arm study this enumerates EVERY unordered pair of arms, not just
  # each arm against the first one. The identified set is spanned by all pairwise
  # arm differences at their shared covariate support, and taking all pairs makes
  # the result independent of the order the arms are listed in.
  if (!is.null(ipd) && nrow(ipd)) {
    for (ss in unique(as.character(ipd[[study]]))) {
      sub <- ipd[as.character(ipd[[study]]) == ss, , drop = FALSE]
      arms <- unique(as.character(sub[[trt]]))
      if (length(arms) < 2L) next
      xof <- function(a) {
        as.matrix(sub[as.character(sub[[trt]]) == a, effect_modifiers,
                      drop = FALSE])
      }
      crow <- function(a) C[match(a, rownames(C)), ]

      # Every PAIR of arms, not just each arm against the study's first arm. An
      # arm contrast a minus b is identified wherever arms a and b BOTH have
      # covariate support, because at a covariate value both attain, the
      # difference of their predictors cancels the study intercept and the
      # prognostic effects and isolates (C_a - C_b)'(beta + Gamma x). Anchoring
      # on the first arm made the verdict depend on which arm a study happened to
      # list first, and it missed contrasts between two arms that share a
      # covariate value when a third arm does not (e.g. a placebo seen only at a
      # different covariate value). Enumerating unordered pairs is
      # order-independent and strictly less conservative, with no false positive:
      # each row is a functional the within-study data genuinely identify.
      for (i in seq_along(arms)) {
        for (j in seq_along(arms)) {
          if (j <= i) next
          m <- crow(arms[i]) - crow(arms[j])
          W <- if (Q) {
            .cpaic_arm_directions(xof(arms[i]), xof(arms[j]), tol = tol)
          } else {
            matrix(1, 1L, 1L)                        # Q = 0: only the level
          }
          for (k in seq_len(ncol(W))) {
            w <- W[, k]
            rows[[length(rows) + 1L]] <-
              as.numeric(kronecker(matrix(w, nrow = 1L), matrix(m, nrow = 1L)))
          }
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
#' @section Strength of the guarantee:
#' The `basis` column states how much the criterion actually proves for each
#' contrast, which is not the same for every row.
#'
#' \describe{
#'   \item{`"exact"`}{The contrast is identified by individual-patient data
#'     under an injective link with no extra support-dependent nuisance, which
#'     means a binomial, poisson or gaussian IPD arm contrast. The IPD likelihood
#'     is then an ordinary regression in arm and covariates, so the within-study
#'     arm-by-covariate variation pins down `m'beta` and `m'Gamma` directly.
#'     Survival is excluded even for IPD, because its flexible baseline hazard
#'     and delayed entry add support-dependent nuisance parameters the
#'     covariate-support argument does not account for.}
#'   \item{`"first-order screen"`}{The contrast is estimable by the linear
#'     row-space criterion, but that criterion is only a design-based screen
#'     here, not an exactness proof, so it can be **optimistic**. This covers two
#'     situations. Identification through aggregate arms under a nonlinear link:
#'     the aggregate likelihood is an integral over the covariate distribution,
#'     and a study pins the contrast down at a variance-weighted mean rather than
#'     at its raw covariate mean. With a log link, one aggregate study and a
#'     symmetric covariate `P(x = -1) = P(x = +1) = 1/2`, the arm means are
#'     `exp(mu)` and `exp(mu + beta) cosh(gamma)`, so the data identify only
#'     `beta + log cosh(gamma)`, not `beta` itself. And identification through
#'     aggregate arms under the identity link: this is exact only when the arms
#'     of each contributing study share a covariate distribution, so that the
#'     study intercept and the prognostic effects cancel from the contrast; cpaic
#'     does not enforce that balance, so it is reported as a screen rather than
#'     claimed exact. Verify these with [prior_sensitivity()].}
#'   \item{`"not identified"`}{Not in the row space of the first-order
#'     information. Any number reported here would be the prior, not the data.
#'     These contrasts are returned as `NA` by [relative_effects()] and dropped
#'     by [cpaic_ranks()].}
#' }
#'
#' @param object A `cpaic_mlnmr` fit.
#' @param newdata A one-row data frame giving the target population's
#'   effect-modifier values. Defaults to the covariate origin.
#' @param reference Reference treatment. Defaults to the fit's reference.
#' @param ... Unused.
#' @return A data frame with `treatment`, `comparator`, `estimable`,
#'   `identified_by` (`"IPD"`, `"aggregate"`, or `"none"`) and `basis`
#'   (`"exact"`, `"first-order screen"`, or `"not identified"`); see the section
#'   below.
#' @references
#' Wigle A, Beliveau A, Nikolakopoulou A, Lin L (2026). Creating Treatment and
#' Component Hierarchies in Component Network Meta-Analysis.
#' @seealso [prior_sensitivity()], [relative_effects()], [cpaic_ranks()]
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

  # "exact" is claimed only where it can be proved, which is narrower than it
  # first appears.
  #
  #  * Individual-patient data under an INJECTIVE link with no extra
  #    support-dependent nuisance: exact. The arm contrast is an ordinary
  #    regression in arm and covariate, so m'beta and m'Gamma are identified
  #    directly. This covers binomial, poisson and gaussian IPD.
  #
  #  * Survival is EXCLUDED even for IPD, because the flexible baseline hazard
  #    and delayed entry add support-dependent nuisance parameters (a study can
  #    fail to identify a baseline interval it has no exposure in) that the
  #    covariate-support argument does not see. Labeling these "exact" overstates
  #    the guarantee.
  #
  #  * Aggregate identification is NOT labeled exact, not even under the identity
  #    link. The identity-link contrast is exact only when the arms of each
  #    contributing aggregate study share a covariate distribution, so that the
  #    study intercept and the prognostic effects cancel. cpaic does not enforce
  #    that balance, and when arm means differ the prognostic coefficient
  #    contaminates the contrast. Reporting such a contrast as a design-based
  #    screen is the honest description.
  #
  # Everything else that is estimable is a first-order screen: correct in rank
  # structure, but able to be optimistic about the anchor. Verify with
  # prior_sensitivity().
  ipd_exact <- by_ipd & !identical(object$family, "survival")
  basis <- ifelse(!ok, "not identified",
                  ifelse(ipd_exact, "exact", "first-order screen"))

  out <- data.frame(
    treatment = others,
    comparator = reference,
    estimable = ok,
    identified_by = ifelse(!ok, "none", ifelse(by_ipd, "IPD", "aggregate")),
    basis = basis,
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(out, "target") <- stats::setNames(x, ems)
  class(out) <- c("cpaic_estimable", "data.frame")
  out
}

#' @export
print.cpaic_estimable <- function(x, ...) {
  tgt <- attr(x, "target")
  cat("Estimability of the population-adjusted relative effects\n")
  if (length(tgt)) {
    cat("  Target population: ",
        paste(names(tgt), signif(tgt, 3), sep = " = ", collapse = ", "),
        "\n", sep = "")
  }
  print(as.data.frame(x), row.names = FALSE)
  if (any(x$basis == "first-order screen")) {
    cat("\n  Rows marked \"first-order screen\" are estimable by the linear",
        " criterion, which\n  is only a design-based screen for them (aggregate",
        " identification, or a\n  survival baseline) and can be optimistic.",
        " Check them with prior_sensitivity().\n", sep = "")
  }
  if (any(x$basis == "not identified")) {
    cat("\n  Rows marked \"not identified\" carry no first-order information;",
        " a number\n  reported for them would be the prior. relative_effects()",
        " returns NA there.\n", sep = "")
  }
  invisible(x)
}
