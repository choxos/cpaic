test_that("fixed design follows column-major gamma ordering", {
  # Given: two components, two covariates, and a nontrivial modifier mapping.
  study <- c(1L, 2L)
  Tc <- matrix(c(1, 2,
                 3, 4), nrow = 2L, byrow = TRUE)
  X <- matrix(c(5, 7,
                11, 13), nrow = 2L, byrow = TRUE)
  gamma <- matrix(c(0.7, 0.8,
                    0.9, 1.0), nrow = 2L)

  # When: the fixed design is built for gamma[, 1] with X[, 2] and
  # gamma[, 2] with X[, 1].
  Z <- cpaic:::.cpaic_fixed_design(
    study = study, Tc = Tc, X = X, n_studies = 2L, emc = c(2L, 1L)
  )

  # Then: each interaction block is component-fast, matching as.vector(gamma).
  interaction <- Z[, 7:10, drop = FALSE]
  expect_equal(
    unname(interaction),
    rbind(c(Tc[1, ] * X[1, 2], Tc[1, ] * X[1, 1]),
          c(Tc[2, ] * X[2, 2], Tc[2, ] * X[2, 1]))
  )

  mu <- c(0.1, 0.2)
  beta <- c(0.3, 0.4)
  breg <- c(0.5, 0.6)
  allbeta <- c(mu, beta, breg, as.vector(gamma))
  expected <- mu[study] + drop(Tc %*% beta) + drop(X %*% breg) +
    drop(Tc %*% gamma[, 1]) * X[, 2] +
    drop(Tc %*% gamma[, 2]) * X[, 1]
  expect_equal(drop(Z %*% allbeta), expected)
})

test_that("thin QR matches multinma scaling and preserves the predictor", {
  # Given: a full column-rank fixed-effects design.
  set.seed(2718)
  Z <- cbind(1, matrix(stats::rnorm(70), nrow = 14L, ncol = 5L))
  allbeta <- seq_len(ncol(Z)) / 10
  qrZ <- qr(Z)
  expected_Q <- qr.Q(qrZ) * sqrt(nrow(Z) - 1)
  expected_R <- qr.R(qrZ)[, sort.list(qrZ$pivot)] /
    sqrt(nrow(Z) - 1)

  # When: the cpaic thin QR helper is applied.
  got <- cpaic:::.cpaic_thin_qr(Z)

  # Then: its factors use the multinma convention and give the same predictor.
  expect_equal(got$Q, expected_Q)
  expect_equal(got$R_inv, solve(expected_R))
  beta_tilde <- solve(got$R_inv, allbeta)
  expect_equal(drop(got$Q %*% beta_tilde), drop(Z %*% allbeta))
})

test_that("QR is an opt-in reparameterization", {
  # Given: the public cmlnmr formals.
  f <- formals(cmlnmr)

  # When: the QR default is inspected.
  default <- eval(f$QR)

  # Then: existing calls retain the original parameterization.
  expect_false(default)
})
