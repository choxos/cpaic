# Shared small anchored data generator: IPD study A vs C with effect
# modification on a binary covariate x, plus an AgD B vs C contrast.
.make_anchored <- function(seed = 42, n = 600, target = 0.30) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.55)
  arm <- rep(c("A", "C"), each = n / 2)
  eta <- -0.2 + 0.8 * (arm == "A") + 0.5 * x - 0.6 * (arm == "A") * x
  y <- rbinom(n, 1, plogis(eta))
  ipd <- data.frame(.study = "S1", .trt = arm, .y = y, x = x)
  agd <- data.frame(studlab = c("S1", "S2"), treat1 = c("A", "B"),
                    treat2 = c("C", "C"), TE = c(NA, 0.4), seTE = c(NA, 0.15))
  list(ipd = ipd, agd = agd, target = target)
}

test_that("cmaic reproduces a manual maicplus anchored Bucher comparison", {
  skip_if_not_installed("maicplus")
  d <- .make_anchored()
  net <- cpaic_network(d$agd, ipd = d$ipd, sm = "OR", family = "binomial",
                       ipd_covariates = "x", reference = "C")
  fit <- cmaic(net, target = c(x = d$target), effect_modifiers = "x",
               n_boot = 50, seed = 1)

  ipd2 <- d$ipd
  ipd2$x_CENTERED <- ipd2$x - d$target
  w <- maicplus::estimate_weights(ipd2, centered_colnames = "x_CENTERED",
                                  boot_strata = ".trt")
  g <- suppressWarnings(stats::glm(.y ~ I(.trt == "A"), family = binomial,
                                   weights = w$data$weights, data = ipd2))
  AvC <- unname(stats::coef(g)[2])

  te <- fit$bridge$fit$TE.random
  expect_equal(unname(te["A", "C"]), AvC, tolerance = 1e-4)
  expect_equal(unname(te["A", "B"]), AvC - 0.4, tolerance = 1e-4)
  expect_true(fit$ess["S1"] > 0 && fit$ess["S1"] <= nrow(d$ipd))
})

test_that("cstc reproduces a manual anchored STC regression", {
  d <- .make_anchored()
  net <- cpaic_network(d$agd, ipd = d$ipd, sm = "OR", family = "binomial",
                       ipd_covariates = "x", reference = "C")
  fit <- suppressWarnings(cstc(net, target = c(x = d$target),
                               effect_modifiers = "x"))

  dd <- d$ipd
  dd$xc <- dd$x - d$target
  dd$.arm <- stats::relevel(factor(dd$.trt), ref = "C")
  g <- stats::glm(.y ~ .arm + xc + .arm:xc, family = binomial, data = dd)
  AvC <- unname(stats::coef(g)[".armA"])

  te <- fit$bridge$fit$TE.random
  expect_equal(unname(te["A", "C"]), AvC, tolerance = 1e-4)
  expect_equal(unname(te["A", "B"]), AvC - 0.4, tolerance = 1e-4)
})

test_that("population adjustment moves the estimate away from naive", {
  d <- .make_anchored()
  dd <- d$ipd
  dd$.arm <- stats::relevel(factor(dd$.trt), ref = "C")
  naive <- unname(stats::coef(stats::glm(.y ~ .arm, family = binomial,
                                         data = dd))[".armA"])
  net <- cpaic_network(d$agd, ipd = d$ipd, sm = "OR", family = "binomial",
                       ipd_covariates = "x", reference = "C")
  fit <- suppressWarnings(cstc(net, target = c(x = d$target),
                               effect_modifiers = "x"))
  adj <- unname(fit$bridge$fit$TE.random["A", "C"])
  expect_false(isTRUE(all.equal(adj, naive, tolerance = 1e-3)))
})
