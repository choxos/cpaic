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
