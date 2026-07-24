# Baseline-shift random-effects correlation.

test_that("the baseline-shift correlation matches multinma::RE_cor", {
  skip_if_not_installed("multinma")

  # Studies listed out of alphabetical order, with a three-arm study, so that a
  # row-ordering mistake in either construction would show up.
  ipd <- data.frame(
    .study = c(rep("z", 6), rep("a", 4)),
    .trt = c(rep(c("C", "A", "B"), each = 2), rep(c("D", "B"), each = 2)),
    stringsAsFactors = FALSE)
  agd <- data.frame(.study = "m", .trt = c("C", "A", "D"),
                    stringsAsFactors = FALSE)

  re <- .cpaic_random_effects(ipd, agd, ".study", ".trt")
  reference <- multinma::RE_cor(
    study = re$arms$study,
    trt = factor(re$arms$trt, levels = sort(unique(re$arms$trt))),
    contrast = rep(FALSE, nrow(re$arms)),
    type = "blshift")

  expect_equal(tcrossprod(re$L_delta), reference, ignore_attr = TRUE)
  expect_equal(nrow(reference), re$N_delta)
})

test_that("the correlation is 1 on the diagonal and 0.5 within a study", {
  R <- .cpaic_re_cor_blshift(c("s1", "s1", "s2"))
  expect_equal(diag(R), rep(1, 3))
  expect_equal(R[1, 2], 0.5)
  expect_equal(R[1, 3], 0)
  # Positive definite, so its Cholesky factor exists.
  expect_silent(chol(R))
})
