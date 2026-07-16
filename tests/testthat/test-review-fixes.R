# Regression tests for the adversarial-review findings.

test_that("interval censoring rejects entry after the interval start", {
  # An interval-censored subject cannot be known to fail after a time it was not
  # yet observed for. entry > start would form (S(start)-S(time))/S(entry) with
  # entry above start, which is not a probability and can exceed one.
  bad <- data.frame(t = 3, s = 3L, st = 1, en = 2)   # entry 2 > start 1
  expect_error(
    cpaic:::.cpaic_survival_outcomes(bad, "t", "s", "st", "en", "test"),
    "entry <= start")
  good <- data.frame(t = 3, s = 3L, st = 2, en = 1)  # entry 1 <= start 2
  expect_silent(
    cpaic:::.cpaic_survival_outcomes(good, "t", "s", "st", "en", "test"))
})

test_that("multi-arm estimability does not depend on arm order", {
  # P is seen only at x = 0; A and B are both seen at x = 1, with enough patients
  # per arm that each arm's covariate support is a genuine (rank-correct) direction
  # rather than an artifact of having fewer rows than parameters. The contrast
  # B vs A at x = 1 is a difference of two observed arm predictors, so it is
  # identified no matter which arm the study lists first; A vs P is not, because
  # those two arms share no covariate value.
  C <- build_C_matrix(c("P", "A", "B"), inactive = "P")
  design <- function(order) {
    ipd <- data.frame(
      .study = "S",
      .trt = rep(c("P", "A", "B"), each = 40),
      x1 = rep(c(0, 1, 1), each = 40))       # constant within arm, 40 each
    ipd <- ipd[order(match(ipd$.trt, order)), ]
    cpaic:::.cpaic_null_space(
      cpaic:::.cpaic_joint_design(C, ipd, NULL, "x1"))
  }
  est <- function(N, t1, t2, x) {
    v <- cpaic:::.cpaic_target_vec(C[t1, ] - C[t2, ], x)
    cpaic:::.cpaic_in_rowspace(matrix(v, nrow = 1L), N)
  }
  for (ord in list(c("P", "A", "B"), c("A", "B", "P"), c("B", "P", "A"))) {
    N <- design(ord)
    expect_true(est(N, "B", "A", 1))     # shared support at x = 1: identified
    expect_false(est(N, "A", "P", 0))    # no shared covariate value: not
  }
})

test_that("a thin arm does not fabricate covariate support", {
  # One patient per arm with two effect modifiers: each arm has fewer rows than
  # the 1 + Q augmented parameters. A recycling bug in the row-space basis used
  # to make such an arm look fully supported, so a prior-only contrast was called
  # identified. With arm perfectly confounded with the covariate, nothing is
  # identifiable, and the joint design must have rank zero.
  C <- build_C_matrix(c("P", "A"), inactive = "P")
  ipd <- data.frame(.study = "S", .trt = c("P", "A"),
                    x1 = c(0, 1), x2 = c(0, 0))
  D <- cpaic:::.cpaic_joint_design(C, ipd, NULL, c("x1", "x2"))
  expect_equal(qr(D)$rank, 0L)
  N <- cpaic:::.cpaic_null_space(D)
  v <- cpaic:::.cpaic_target_vec(C["A", ] - C["P", ], c(0, 0))
  expect_false(cpaic:::.cpaic_in_rowspace(matrix(v, nrow = 1L), N))
})

test_that(".cpaic_row_space returns the true rank for a wide matrix", {
  # A single augmented row (1, 0) spans a one-dimensional space, not two.
  expect_equal(ncol(cpaic:::.cpaic_row_space(cbind(1, 0))), 1L)
  # Two identical rows still span one dimension.
  expect_equal(ncol(cpaic:::.cpaic_row_space(rbind(c(1, 2, 3), c(1, 2, 3)))), 1L)
})

test_that("edge_influence tolerance is relative to the largest influence", {
  # Three edges of equal influence 1/3. At tol = 0.5 the documented cutoff is
  # 0.5 * (1/3), so none is "dead". An absolute floor of 1 in the scale would
  # make the cutoff 0.5 and fire a false warning on every edge.
  net <- cpaic_network(
    data.frame(
      studlab = c("s1", "s2", "s3"),
      treat1 = c("A", "A", "A"),
      treat2 = c("B", "B", "B"),
      TE = c(0.5, 0.5, 0.5), seTE = c(0.2, 0.2, 0.2)),
    sm = "OR", inactive = "A")
  br <- cnma_bridge(net)
  expect_silent(ei <- edge_influence(br, treatment = "B", tol = 0.5))
  expect_true(all(abs(ei$influence) > 0))
})

test_that("copula correlation is adjusted to the latent scale for binaries", {
  # multinma's Pearson adjustment: sin(pi r / 2) for binary-binary, sqrt(pi/2) r
  # for binary-continuous, continuous-continuous unchanged.
  bb <- cpaic:::.cpaic_cor_adjust_pearson(
    matrix(c(1, 0.5, 0.5, 1), 2), c("bernoulli", "bernoulli"))
  expect_equal(bb[1, 2], sin(pi * 0.5 / 2))
  cb <- cpaic:::.cpaic_cor_adjust_pearson(
    matrix(c(1, 0.4, 0.4, 1), 2), c("normal", "bernoulli"))
  expect_equal(cb[1, 2], sqrt(pi / 2) * 0.4)
  nn <- cpaic:::.cpaic_cor_adjust_pearson(
    matrix(c(1, 0.6, 0.6, 1), 2), c("normal", "normal"))
  expect_equal(nn[1, 2], 0.6)
})
