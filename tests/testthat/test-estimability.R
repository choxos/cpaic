# Estimability of population-adjusted component contrasts.
#
# The criterion extends the row-space result of Wigle et al. (2026) from beta to
# (beta, vec(Gamma)). These tests pin down the properties it must have, including
# two counterexamples that an earlier version got wrong.

C3 <- build_C_matrix(c("Placebo", "A", "B", "A+B", "A+B+C"),
                     inactive = "Placebo")

est_at <- function(D, m, x) {
  N <- cpaic:::.cpaic_null_space(D)
  cpaic:::.cpaic_in_rowspace(
    matrix(cpaic:::.cpaic_target_vec(m, x), nrow = 1L), N)
}

test_that("Q = 0 collapses to the Wigle et al. (2026) row-space criterion", {
  Cq <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")

  # (a) Both edges observed: A and A+B vs Placebo. Then B = (A+B) - A is
  #     identified as the difference, so every component effect is estimable.
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 10))
  agd <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"))
  D <- cpaic:::.cpaic_joint_design(Cq, ipd, agd, character(0))
  expect_equal(ncol(D), ncol(Cq))              # Q = 0 -> phi = beta
  expect_true(est_at(D, Cq["A+B", ] - Cq["Placebo", ], numeric(0)))
  expect_true(est_at(D, Cq["A", ] - Cq["Placebo", ], numeric(0)))
  expect_true(est_at(D, Cq["A+B", ] - Cq["A", ], numeric(0)))   # = B

  # (b) Only A+B vs Placebo observed. The SUM is identified, but neither
  #     component separately: the canonical Wigle et al. (2026) example.
  agd2 <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"))
  D2 <- cpaic:::.cpaic_joint_design(Cq, NULL, agd2, character(0))
  expect_true(est_at(D2, Cq["A+B", ] - Cq["Placebo", ], numeric(0)))  # the sum
  expect_false(est_at(D2, Cq["A", ] - Cq["Placebo", ], numeric(0)))   # not A
  expect_false(est_at(D2, Cq["A+B", ] - Cq["A", ], numeric(0)))       # not B
})

test_that("an IPD study with a CONSTANT modifier is anchored at its own value", {
  # Grok counterexample A. If x1 never varies within the trial (a single-sex
  # trial, say) then only m'(beta + 5 Gamma) is learned. Anchoring the level row
  # at the covariate origin would claim x1 = 0 is estimable (it is not) and deny
  # x1 = 5 (which is).
  Cc <- build_C_matrix(c("Placebo", "A"), inactive = "Placebo")
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 20),
                    x1 = 5)
  D <- cpaic:::.cpaic_joint_design(Cc, ipd, NULL, "x1")
  m <- Cc["A", ] - Cc["Placebo", ]

  expect_equal(nrow(D), 1L)          # level only; no slope row
  expect_true(est_at(D, m, 5))       # estimable exactly where it was held
  expect_false(est_at(D, m, 0))      # NOT at the covariate origin
  expect_false(est_at(D, m, 1))
})

test_that("partial covariate variation transports only in the varying direction", {
  # x1 varies, x2 is held at 5. We may transport in x1 but never in x2.
  Cc <- build_C_matrix(c("Placebo", "A"), inactive = "Placebo")
  set.seed(4)
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 30),
                    x1 = rnorm(60), x2 = 5)
  D <- cpaic:::.cpaic_joint_design(Cc, ipd, NULL, c("x1", "x2"))
  m <- Cc["A", ] - Cc["Placebo", ]

  expect_true(est_at(D, m, c(0, 5)))     # any x1, but x2 must be 5
  expect_true(est_at(D, m, c(2, 5)))
  expect_false(est_at(D, m, c(0, 0)))    # x2 = 0 was never observed
  expect_false(est_at(D, m, c(0, 1)))
})

test_that("an IPD study with full covariate variation transports anywhere", {
  set.seed(5)
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 40),
                    x1 = rnorm(80))
  D <- cpaic:::.cpaic_joint_design(C3, ipd, NULL, "x1")
  m <- C3["A", ] - C3["Placebo", ]
  for (x in c(-2, 0, 0.5, 3)) expect_true(est_at(D, m, x))
})

test_that("an aggregate contrast is estimable only in its own population", {
  # One aggregate two-arm study gives ONE equation for 1 + Q unknowns, so it
  # pins the contrast down at its own covariate mean and nowhere else.
  agd <- data.frame(.study = "S2", .trt = c("Placebo", "B"),
                    x1_mean = c(0.5, 0.5))
  D <- cpaic:::.cpaic_joint_design(C3, NULL, agd, "x1")
  m <- C3["B", ] - C3["Placebo", ]

  expect_true(est_at(D, m, 0.5))      # its own mean
  expect_false(est_at(D, m, 0))       # anywhere else
  expect_false(est_at(D, m, 1))
})

test_that("two aggregate studies span the line through their means", {
  # Ecological identification: differing study means give a between-study
  # gradient that spans the interaction direction, so the whole affine hull of
  # the means becomes estimable.
  agd <- data.frame(.study = c("S2", "S2", "S3", "S3"),
                    .trt = c("Placebo", "B", "Placebo", "B"),
                    x1_mean = c(0.5, 0.5, -0.5, -0.5))
  D <- cpaic:::.cpaic_joint_design(C3, NULL, agd, "x1")
  m <- C3["B", ] - C3["Placebo", ]
  for (x in c(-0.5, 0, 0.5, 2)) expect_true(est_at(D, m, x))
})

test_that("population adjustment is strictly harder than reconnection", {
  # The cross-gap contrast A+B+C vs Placebo IS estimable as an aggregate-data
  # component contrast, but is NOT estimable as a population-adjusted contrast,
  # because Gamma_B and Gamma_C are never identified.
  set.seed(6)
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 40),
                    x1 = rnorm(80))
  agd <- data.frame(.study = c("S2", "S2", "S3", "S3"),
                    .trt = c("Placebo", "B", "A+B", "A+B+C"),
                    x1_mean = c(0.5, 0.5, -0.4, -0.4))
  D <- cpaic:::.cpaic_joint_design(C3, ipd, agd, "x1")
  m <- C3["A+B+C", ] - C3["Placebo", ]

  expect_false(est_at(D, m, 1))
  expect_false(est_at(D, m, 0))

  # The Q = 0 (aggregate cNMA) view says it IS estimable, which is exactly the
  # trap: reconnecting the network does not license population adjustment.
  D0 <- cpaic:::.cpaic_joint_design(C3, ipd, agd, character(0))
  expect_true(est_at(D0, m, numeric(0)))
})

test_that("multi-arm IPD contributes its full contrast span", {
  # All arm contrasts against the first arm are included, so the span is the
  # study's whole contrast space. Keeping only one edge would understate it.
  set.seed(7)
  ipd <- data.frame(.study = "S1",
                    .trt = rep(c("Placebo", "A", "B"), each = 30),
                    x1 = rnorm(90))
  D <- cpaic:::.cpaic_joint_design(C3, ipd, NULL, "x1")
  # A vs B is estimable even though neither is the study's reference arm.
  expect_true(est_at(D, C3["A", ] - C3["B", ], 1.5))
})
