# Fixed-effects design and thin QR reparameterization for cML-NMR.
# The QR scaling and pivot restoration follow multinma, which implements the
# ML-NMR methods of Phillippo et al. (2020). Both packages are GPL-3.

#' Build the complete cML-NMR fixed-effects design
#'
#' The interaction columns are ordered by effect modifier, then component.
#' This is the column-major order used by `as.vector(gamma)` in R and
#' `to_vector(gamma)` in Stan.
#' @noRd
.cpaic_fixed_design <- function(study, Tc, X, n_studies,
                                emc = seq_len(ncol(X))) {
  Tc <- as.matrix(Tc)
  X <- as.matrix(X)
  n <- nrow(Tc)
  if (nrow(X) != n || length(study) != n) {
    stop("`study`, `Tc`, and `X` must describe the same number of rows.",
         call. = FALSE)
  }
  if (length(emc) == 0L || any(emc < 1L) || any(emc > ncol(X))) {
    stop("`emc` must select valid columns of `X`.", call. = FALSE)
  }
  if (any(study < 1L) || any(study > n_studies)) {
    stop("`study` contains an index outside `n_studies`.", call. = FALSE)
  }

  study_design <- matrix(0, nrow = n, ncol = n_studies)
  study_design[cbind(seq_len(n), study)] <- 1
  interaction <- do.call(
    cbind,
    lapply(emc, function(q) sweep(Tc, 1L, X[, q], `*`))
  )

  colnames(study_design) <- paste0("mu[", seq_len(n_studies), "]")
  colnames(Tc) <- paste0("beta[", seq_len(ncol(Tc)), "]")
  colnames(X) <- paste0("breg[", seq_len(ncol(X)), "]")
  colnames(interaction) <- unlist(
    lapply(seq_along(emc), function(q) {
      paste0("gamma[", seq_len(ncol(Tc)), ",", q, "]")
    }),
    use.names = FALSE
  )

  cbind(study_design, Tc, X, interaction)
}

#' Compute the scaled thin QR decomposition used by multinma
#' @noRd
.cpaic_thin_qr <- function(Z) {
  Z <- as.matrix(Z)
  if (nrow(Z) <= 1L || nrow(Z) < ncol(Z)) {
    stop("QR requires at least two rows and no fewer rows than columns.",
         call. = FALSE)
  }

  qrZ <- qr(Z)
  if (qrZ$rank < ncol(Z)) {
    stop("QR requires a full column-rank fixed-effects design.",
         call. = FALSE)
  }
  Q <- qr.Q(qrZ) * sqrt(nrow(Z) - 1)
  R <- qr.R(qrZ)[, sort.list(qrZ$pivot), drop = FALSE] /
    sqrt(nrow(Z) - 1)

  list(Q = Q, R_inv = solve(R))
}
