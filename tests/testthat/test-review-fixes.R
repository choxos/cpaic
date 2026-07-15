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
  # P is seen only at x = 0; A and B are both seen at x = 1. The contrast B vs A
  # at x = 1 is a difference of two observed arm predictors, so it is identified
  # no matter which arm the study happens to list first.
  C <- build_C_matrix(c("P", "A", "B"), inactive = "P")
  ident_B_vs_A <- function(order) {
    ipd <- data.frame(
      .study = "S",
      .trt = c("P", "A", "B"),
      x1 = c(0, 1, 1))
    ipd <- ipd[order(match(ipd$.trt, order)), ]
    N <- cpaic:::.cpaic_null_space(
      cpaic:::.cpaic_joint_design(C, ipd, NULL, "x1"))
    v <- cpaic:::.cpaic_target_vec(C["B", ] - C["A", ], 1)
    cpaic:::.cpaic_in_rowspace(matrix(v, nrow = 1L), N)
  }
  expect_true(ident_B_vs_A(c("P", "A", "B")))
  expect_true(ident_B_vs_A(c("A", "B", "P")))
  expect_true(ident_B_vs_A(c("B", "P", "A")))
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
