# Regression tests for correctness fixes and identifiability edge cases.

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

# --- Data-integrity and fit-validity guards ----------------------------------

test_that("replace_contrasts rejects duplicate keys and does not clone metadata", {
  cols <- list(studlab = "studlab", treat1 = "treat1", treat2 = "treat2",
               TE = "TE", seTE = "seTE")
  # Two aggregate rows for the same study and unordered pair would each be
  # overwritten with the same adjusted estimate and double-counted.
  dup <- data.frame(studlab = c("s1", "s1"), treat1 = c("A", "A"),
                    treat2 = c("B", "B"), TE = c(0.1, 0.1), seTE = c(0.2, 0.2),
                    extra = c("x", "y"), stringsAsFactors = FALSE)
  adj <- data.frame(studlab = "s1", treat1 = "A", treat2 = "B",
                    TE = 0.5, seTE = 0.1, stringsAsFactors = FALSE)
  expect_error(cpaic:::.cpaic_replace_contrasts(dup, adj, cols), "duplicate")

  # An appended edge must not inherit unrelated metadata from row 1.
  agd <- data.frame(studlab = "s2", treat1 = "C", treat2 = "D",
                    TE = 0.3, seTE = 0.2, extra = "keepme",
                    stringsAsFactors = FALSE)
  out <- cpaic:::.cpaic_replace_contrasts(agd, adj, cols)
  newrow <- out[out$studlab == "s1", ]
  expect_true(is.na(newrow$extra))            # not cloned from agd[1, ]
  expect_equal(newrow$TE, 0.5)
  expect_equal(out$extra[out$studlab == "s2"], "keepme")
})

test_that("a named copula correlation is reordered to the effect modifiers", {
  em <- c("x1", "x2", "x3")
  R <- matrix(c(1, 0.2, 0.4, 0.2, 1, 0.6, 0.4, 0.6, 1), 3)
  dimnames(R) <- list(c("x3", "x1", "x2"), c("x3", "x1", "x2"))  # scrambled
  out <- cpaic:::.cpaic_copula_cor(NULL, em, ".study", given = R)
  # In input order (x3, x1, x2): cor(x1, x2) = R[2,3] = 0.6, cor(x1, x3) = 0.2,
  # cor(x2, x3) = 0.4. Reordered to (x1, x2, x3):
  expect_equal(out[1, 2], 0.6)
  expect_equal(out[1, 3], 0.2)
  expect_equal(out[2, 3], 0.4)
  Rbad <- R
  dimnames(Rbad) <- list(c("z", "x1", "x2"), c("z", "x1", "x2"))
  expect_error(cpaic:::.cpaic_copula_cor(NULL, em, ".study", given = Rbad),
               "dimnames")
})

test_that("construction rejects self-comparisons, duplicates, and bad coding", {
  expect_error(
    cpaic_network(data.frame(studlab = "s1", treat1 = "A", treat2 = "A",
                             TE = 0.1, seTE = 0.2), sm = "OR"),
    "self-comparison")
  expect_error(
    cpaic_network(data.frame(studlab = c("s1", "s1"), treat1 = c("A", "B"),
                             treat2 = c("B", "A"), TE = c(0.1, 0.1),
                             seTE = c(0.2, 0.2)), sm = "OR"),
    "duplicate")
  expect_error(build_C_matrix(c("A", "A++B")), "empty component")
  expect_error(build_C_matrix(c("A", "B"), inactive = "Nope"), "no-op")
})

test_that("the bridge rejects the empty common = random = FALSE model", {
  net <- cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
  expect_error(cnma_bridge(net, common = FALSE, random = FALSE),
               "At least one")
})

test_that("n_boot must be a whole number", {
  net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                       family = "binomial", ipd_covariates = "x1",
                       inactive = "Placebo")
  expect_error(cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
                     n_boot = 2.5), ">= 2")
})

test_that("a separated two-stage fit is rejected, not passed to the bridge", {
  # Complete separation: every Placebo patient a 0, every A patient a 1.
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 20),
                    .y = c(rep(0, 20), rep(1, 20)), x1 = stats::rnorm(40))
  agd <- data.frame(studlab = c("S1", "S2"), treat1 = c("A", "A+B"),
                    treat2 = c("Placebo", "Placebo"), TE = c(0.5, 0.6),
                    seTE = c(0.2, 0.2))
  net <- cpaic_network(agd, ipd = ipd, sm = "OR", family = "binomial",
                       ipd_covariates = "x1", inactive = "Placebo")
  expect_error(cstc(net, target = c(x1 = 0), effect_modifiers = "x1"),
               "not usable|separation")
})

test_that("the estimability screen is invariant to aggregate row multiplicity", {
  C <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")
  agd1 <- data.frame(.study = "S", .trt = c("Placebo", "A"),
                     x1_mean = c(0.2, 0.2))
  # The same two arms, but each arm's summary repeated a different number of
  # times (as survival pseudo-IPD would store it).
  agd_rep <- data.frame(
    .study = "S",
    .trt = c(rep("Placebo", 5), rep("A", 50)),
    x1_mean = 0.2)
  D1 <- cpaic:::.cpaic_joint_design(C, NULL, agd1, "x1")
  D2 <- cpaic:::.cpaic_joint_design(C, NULL, agd_rep, "x1")
  expect_equal(D1, D2)
})

test_that("cmlnmr rejects data overrides, source overlap, and non-integer n_int", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("randtoolbox")
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 10),
                    .y = stats::rbinom(20, 1, 0.5), x1 = stats::rnorm(20))
  agd_ok <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"),
                       r = c(5, 6), n = c(10, 10),
                       x1_mean = c(0, 0), x1_sd = c(1, 1))
  # Overlapping study id in both sources.
  agd_dup <- data.frame(.study = "S1", .trt = c("Placebo", "A+B"),
                        r = c(5, 6), n = c(10, 10),
                        x1_mean = c(0, 0), x1_sd = c(1, 1))
  expect_error(cmlnmr(ipd, agd_dup, effect_modifiers = "x1",
                      inactive = "Placebo"), "both")
  expect_error(cmlnmr(ipd, agd_ok, effect_modifiers = "x1",
                      inactive = "Placebo", data = list(x = 1)),
               "cannot be passed")
  expect_error(cmlnmr(ipd, agd_ok, effect_modifiers = "x1",
                      inactive = "Placebo", n_int = 64.5),
               "positive integer")
  expect_error(cmlnmr(ipd, agd_ok, effect_modifiers = "x1",
                      inactive = "Placebo", seed = 1e12),
               "whole number")
})

test_that("redact_fit removes raw data and blocks refitting", {
  obj <- structure(list(refit_args = list(ipd = data.frame(x = 1),
                                          agd = data.frame(y = 1)),
                        observed = list(1:3)),
                   class = c("cpaic_mlnmr", "cpaic_fit"))
  r <- redact_fit(obj)
  expect_null(r$refit_args$ipd)
  expect_null(r$observed)
  expect_true(isTRUE(attr(r, "redacted")))
  expect_error(prior_sensitivity(r, newdata = data.frame(x = 0)), "redacted")
})

test_that("the fit gate is specific: valid strong effects are not rejected", {
  set.seed(1)
  # Gaussian, large outcome scale, strong real arm effect: not a separation, so
  # the magnitude heuristic must not fire on the unbounded mean-difference scale.
  dg <- data.frame(y = c(stats::rnorm(60, 5000, 300),
                         stats::rnorm(60, 7000, 300)),
                   arm = relevel(factor(rep(c("P", "A"), each = 60)), "P"))
  fg <- stats::glm(y ~ arm, stats::gaussian(), dg)
  expect_length(cpaic:::.cpaic_regression_problems(fg, "gaussian", "armA"), 0L)

  # Binomial with a strong prognostic pushing some fitted values to the 0/1
  # boundary, but a moderate, well-identified treatment effect: not separation.
  x <- c(stats::rnorm(80, -3, 1), stats::rnorm(80, 3, 1))
  y <- stats::rbinom(160, 1, stats::plogis(2 * x))
  arm <- relevel(factor(rep(c("P", "A"), each = 80)), "P")
  fb <- suppressWarnings(stats::glm(y ~ arm + x, stats::binomial()))
  expect_length(cpaic:::.cpaic_regression_problems(fb, "binomial", "armA"), 0L)

  # Complete separation is still caught.
  ys <- c(rep(0, 40), rep(1, 40))
  arm2 <- relevel(factor(rep(c("P", "A"), each = 40)), "P")
  fs <- suppressWarnings(stats::glm(ys ~ arm2, stats::binomial()))
  expect_gt(length(cpaic:::.cpaic_regression_problems(fs, "binomial", "arm2A")),
            0L)
})
