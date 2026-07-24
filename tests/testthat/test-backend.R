# The two sampler backends must be interchangeable.

test_that("backend selection is validated and reported", {
  expect_error(
    cmlnmr(data.frame(), data.frame(), effect_modifiers = "x",
           backend = "nope"),
    "should be one of")
})

test_that("a fit reports which engine produced it", {
  # Dispatch is on the fit's class, not on a recorded flag, so a saved object
  # still resolves correctly.
  expect_identical(.cpaic_fit_backend(structure(list(), class = "stanfit")),
                   "rstan")
  expect_identical(.cpaic_fit_backend(structure(list(), class = "CmdStanMCMC")),
                   "cmdstanr")
})

test_that("cmlnmr fits with the rstan backend", {
  skip_on_cran()

  set.seed(1)
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 40),
                    .y = rbinom(80, 1, 0.5), x1 = rnorm(80))
  agd <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"),
                    r = c(18, 26), n = c(40, 40),
                    x1_mean = c(0.2, 0.2), x1_sd = c(1, 1))
  fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
                backend = "rstan", chains = 2, iter_warmup = 300,
                iter_sampling = 300, seed = 1)

  expect_s4_class(fit$fit, "stanfit")
  expect_identical(fit$backend, "rstan")
  expect_identical(fit$provenance$backend, "rstan")
  # A appears in the IPD contrast; B only inside the aggregate A+B arm, so at
  # the covariate origin the design identifies the sum and not B alone.
  expect_identical(fit$components$component, c("A", "B"))
  expect_true(is.finite(fit$components$estimate[1]))
  expect_true(is.na(fit$components$estimate[2]))
  # The diagnostics must come back as numbers, not the infinities that
  # max()/min() return over an all-NA column.
  expect_true(is.finite(fit$diagnostics$max_rhat))
  expect_true(is.finite(fit$diagnostics$min_ess))

  re <- relative_effects(fit, newdata = data.frame(x1 = 0.2))
  expect_true("basis" %in% names(re))
  expect_true(all(is.finite(re$estimate)))
})

test_that("the two backends agree up to Monte Carlo error", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "cmdstan not installed")

  # rstan and cmdstanr do not share a random number stream, so the same seed
  # gives different draws. What must agree is the posterior they target, and
  # everything the package derives from a fit.
  set.seed(1)
  ipd <- data.frame(.study = "S1", .trt = rep(c("Placebo", "A"), each = 200),
                    .y = rbinom(400, 1, 0.5), x1 = rnorm(400))
  agd <- data.frame(.study = "S2", .trt = c("Placebo", "A+B"),
                    r = c(90, 130), n = c(200, 200),
                    x1_mean = c(0.2, 0.2), x1_sd = c(1, 1))
  args <- list(ipd = ipd, agd = agd, effect_modifiers = "x1",
               inactive = "Placebo", chains = 2, iter_warmup = 600,
               iter_sampling = 600, seed = 1)
  f_r <- do.call(cmlnmr, c(args, list(backend = "rstan")))
  f_c <- do.call(cmlnmr, c(args, list(backend = "cmdstanr")))

  nd <- data.frame(x1 = 0.2)
  r_r <- relative_effects(f_r, newdata = nd)
  r_c <- relative_effects(f_c, newdata = nd)

  expect_identical(r_r$treatment, r_c$treatment)
  expect_identical(r_r$basis, r_c$basis)
  expect_identical(is.na(f_r$components$estimate),
                   is.na(f_c$components$estimate))
  # Generous but not vacuous: these are separate short chains, so agreement
  # within a fraction of the posterior standard deviation is the claim.
  expect_equal(r_r$estimate, r_c$estimate, tolerance = 0.15)
})
