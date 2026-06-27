# Input validation and robustness guards.

test_that("family x sm mismatch is rejected at network construction", {
  expect_error(
    cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "RR",
                  family = "binomial", ipd_covariates = "x1",
                  inactive = "Placebo"),
    "must be one of"
  )
})

test_that("missing IPD outcome column is caught", {
  ipd2 <- cpaic_bin_ipd
  ipd2$.y <- NULL
  expect_error(
    cpaic_network(cpaic_bin_agd, ipd = ipd2, sm = "OR", family = "binomial",
                  ipd_covariates = "x1", inactive = "Placebo"),
    "outcome column"
  )
})

test_that("non-finite target is rejected", {
  net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                       family = "binomial", ipd_covariates = "x1",
                       inactive = "Placebo")
  expect_error(cstc(net, target = c(x1 = NA), effect_modifiers = "x1"),
               "finite numeric")
  expect_error(cmaic(net, target = c(x1 = Inf), effect_modifiers = "x1"),
               "finite numeric")
})

test_that("multi-arm IPD studies are rejected", {
  set.seed(1)
  ipd3 <- data.frame(.study = "S", .trt = rep(c("A", "B", "C"), each = 20),
                     .y = rbinom(60, 1, 0.5), x1 = rnorm(60))
  agd3 <- data.frame(studlab = c("S", "S"), treat1 = c("A", "B"),
                     treat2 = c("C", "C"), TE = c(0.2, 0.1),
                     seTE = c(0.1, 0.1))
  net3 <- cpaic_network(agd3, ipd = ipd3, sm = "OR", family = "binomial",
                        ipd_covariates = "x1")
  expect_error(cstc(net3, target = c(x1 = 0), effect_modifiers = "x1"),
               "two-arm")
})

test_that("cnma_bridge warns on a non-identifiable component network", {
  agdni <- data.frame(studlab = "S1", treat1 = "A+B", treat2 = "Placebo",
                      TE = 0.5, seTE = 0.1)
  netni <- cpaic_network(agdni, sm = "OR", inactive = "Placebo")
  expect_false(cpaic_connectivity(netni)$identifiable)
  expect_warning(cnma_bridge(netni), "not uniquely identifiable")
})

test_that("disconnected + non-identifiable network errors in cnma_bridge", {
  agd <- data.frame(studlab = c("S1", "S2"), treat1 = c("A", "C"),
                    treat2 = c("PlaceboX", "PlaceboY"), TE = c(0.5, 0.4),
                    seTE = c(0.1, 0.1))
  net <- cpaic_network(agd, sm = "OR")
  expect_false(cpaic_connectivity(net)$identifiable)
  expect_error(cnma_bridge(net), "disconnected and cannot be bridged")
})

test_that("empty agd and out-of-network IPD arms are rejected", {
  expect_error(cpaic_network(cpaic_bin_agd[0, ], sm = "OR"), "no comparisons")
  set.seed(1)
  ipd_bad <- data.frame(.study = "S1", .trt = rep(c("Z", "Placebo"), each = 10),
                        .y = rbinom(20, 1, 0.5), x1 = rnorm(20))
  expect_error(
    cpaic_network(cpaic_bin_agd, ipd = ipd_bad, sm = "OR",
                  family = "binomial", ipd_covariates = "x1",
                  inactive = "Placebo"),
    "not present in the aggregate network")
})

test_that("cmaic validates n_boot and target_sd", {
  net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                       family = "binomial", ipd_covariates = "x1",
                       inactive = "Placebo")
  expect_error(cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
                     n_boot = 0), "positive integer")
  expect_error(cmaic(net, target = c(x1 = 0), target_sd = c(x1 = -1),
                     effect_modifiers = "x1", n_boot = 10), "non-negative")
})

test_that("relative_effects validates reference and level", {
  br <- cnma_bridge(cpaic_network(cpaic_bin_agd, sm = "OR",
                                  inactive = "Placebo"))
  expect_error(relative_effects(br, reference = "NOPE"), "must be one of")
  expect_error(relative_effects(br, level = 1.5), "in \\(0, 1\\)")
})

test_that("target_sd matches second moments (no longer a no-op)", {
  net <- cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                       family = "binomial", ipd_covariates = "x1",
                       inactive = "Placebo")
  f1 <- cmaic(net, target = c(x1 = 0), effect_modifiers = "x1",
              n_boot = 10, seed = 1)
  f2 <- cmaic(net, target = c(x1 = 0), target_sd = c(x1 = 0.5),
              effect_modifiers = "x1", n_boot = 10, seed = 1)
  expect_false(isTRUE(all.equal(unname(f1$ess), unname(f2$ess))))
})
