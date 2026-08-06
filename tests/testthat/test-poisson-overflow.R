make_poisson_draw_fit <- function(draw_map) {
  list(draws = function(variables, format) {
    value <- draw_map[[variables]]
    if (is.null(value)) stop("unknown draw variable: ", variables)
    if (is.null(dim(value))) value <- matrix(value, ncol = 1L)
    value
  })
}

make_prior_poisson_fit <- function(draw_map, overflow_variables = TRUE) {
  out <- list(
    fit = make_poisson_draw_fit(draw_map),
    family = "poisson",
    prior_predictive = TRUE,
    observed = list(ipd = c(0L, 2L), agd = 3L),
    rep_variables = list(ipd = "yrep_ipd", agd = "rrep_agd")
  )
  if (overflow_variables) {
    out$rep_overflow_variables <- list(
      ipd = "yrep_ipd_overflow_count",
      agd = "rrep_agd_overflow_count"
    )
  }
  structure(out, class = c("cpaic_mlnmr", "cpaic_fit"))
}

test_that("Poisson fitted means average rates on the response scale", {
  draw_map <- list(
    eta_ipd = matrix(
      c(log(1), log(1), log(9), log(9)), nrow = 2L, byrow = TRUE
    ),
    log_lambda_agd = matrix(c(log(1), log(9)), ncol = 1L)
  )
  object <- structure(
    list(
      fit = make_poisson_draw_fit(draw_map),
      family = "poisson",
      refit_args = list(
        ipd = data.frame(.y = c(0L, 1L)),
        agd = data.frame(r = 2L),
        outcome = ".y",
        r = "r"
      )
    ),
    class = c("cpaic_mlnmr", "cpaic_fit")
  )

  fitted <- cpaic:::.cpaic_fitted_observed(object)$fitted
  expect_equal(fitted, c(5, 5, 5), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(fitted, c(3, 3, 3))))
})

test_that("a real Poisson fit records safe replicated draws without caps", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  skip_if_not_installed("randtoolbox")

  Cmat <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")
  set.seed(71)
  make_arm <- function(study, treatment, n) {
    x1 <- stats::rnorm(n)
    component_effect <- sum(Cmat[treatment, ] * c(A = 0.3, B = 0.2))
    rate <- exp(-1 + 0.15 * x1 + component_effect)
    data.frame(
      .study = study, .trt = treatment,
      .y = stats::rpois(n, rate), .exposure = 1, x1 = x1
    )
  }
  ipd <- rbind(
    make_arm("S1", "Placebo", 40),
    make_arm("S1", "A", 40)
  )
  agd <- data.frame(
    .study = "S2", .trt = c("A", "A+B"),
    r = c(18L, 24L), E = c(40, 40),
    x1_mean = c(0, 0), x1_sd = c(1, 1)
  )

  fit <- suppressWarnings(cmlnmr(
    ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
    family = "poisson", chains = 1, iter_warmup = 100,
    iter_sampling = 100, n_int = 8, seed = 19
  ))

  expect_false(fit$replication_overflow$any)
  expect_true(all(fit$replication_overflow$summary$overflowed_values == 0L))
  expect_true(all(cpaic:::.cpaic_draws(fit, "yrep_ipd") >= 0L))
  expect_true(all(cpaic:::.cpaic_draws(fit, "rrep_agd") >= 0L))
})

test_that("legacy Poisson fits fail closed before predictive summaries", {
  object <- make_prior_poisson_fit(
    list(
      yrep_ipd = matrix(c(0L, 1L, 1L, 2L), nrow = 2L),
      rrep_agd = matrix(c(2L, 3L), ncol = 1L)
    ),
    overflow_variables = FALSE
  )

  expect_error(
    prior_predictive_check(object),
    "overflow metadata are missing or malformed"
  )
})

test_that("Poisson predictive checks reject recorded RNG overflow", {
  object <- make_prior_poisson_fit(list(
    yrep_ipd = matrix(c(0L, 1L, 1L, -1L), nrow = 2L),
    rrep_agd = matrix(c(2L, 3L), ncol = 1L),
    yrep_ipd_overflow_count = c(0L, 1L),
    rrep_agd_overflow_count = c(0L, 0L)
  ))

  expect_error(
    prior_predictive_check(object),
    "RNG overflow sentinels.*ipd=1"
  )
})

test_that("Poisson predictive checks verify sentinels against zero counts", {
  object <- make_prior_poisson_fit(list(
    yrep_ipd = matrix(c(0L, 1L, 1L, -1L), nrow = 2L),
    rrep_agd = matrix(c(2L, 3L), ncol = 1L),
    yrep_ipd_overflow_count = c(0L, 0L),
    rrep_agd_overflow_count = c(0L, 0L)
  ))

  expect_error(
    prior_predictive_check(object),
    "negative overflow sentinels"
  )
})

test_that("valid Poisson replicated outcomes retain predictive summaries", {
  object <- make_prior_poisson_fit(list(
    yrep_ipd = matrix(c(0L, 1L, 1L, 2L), nrow = 2L, byrow = TRUE),
    rrep_agd = matrix(c(2L, 4L), ncol = 1L),
    yrep_ipd_overflow_count = c(0L, 0L),
    rrep_agd_overflow_count = c(0L, 0L)
  ))

  out <- prior_predictive_check(object)
  expect_s3_class(out, "cpaic_prior_predictive")
  expect_equal(out$source, c("ipd", "agd"))
  expect_true(all(is.finite(out$rep_median)))
})

test_that("partial sampler diagnostics remain explicitly unavailable", {
  testthat::local_mocked_bindings(
    .cpaic_sampler_diagnostics = function(fit) list(
      num_divergent = NA_integer_,
      num_max_treedepth = 0L,
      ebfmi = 0.8,
      unavailable = FALSE
    ),
    .cpaic_stan_variables = function(fit) "beta",
    .cpaic_fit_summary = function(fit, variables) data.frame(
      rhat = 1, ess_bulk = 200, ess_tail = 200
    ),
    .package = "cpaic"
  )

  expect_warning(
    out <- cpaic:::.cpaic_check_diagnostics(list()),
    "diagnostics unavailable"
  )
  expect_true(out$unavailable)
  expect_true(is.na(out$divergences))
  expect_equal(out$max_treedepth, 0L)
  expect_identical(out$status, "unknown")
  expect_false(out$decision_grade)
})

test_that("sampler diagnostic validity is a first-class fit state", {
  testthat::local_mocked_bindings(
    .cpaic_sampler_diagnostics = function(fit) list(
      num_divergent = 0L, num_max_treedepth = 0L, ebfmi = 0.8,
      unavailable = FALSE
    ),
    .cpaic_stan_variables = function(fit) "beta",
    .cpaic_fit_summary = function(fit, variables) data.frame(
      rhat = 1.001, ess_bulk = 500, ess_tail = 450
    ),
    .package = "cpaic"
  )
  passed <- cpaic:::.cpaic_check_diagnostics(list())
  expect_identical(passed$status, "passed")
  expect_true(passed$decision_grade)

  testthat::local_mocked_bindings(
    .cpaic_sampler_diagnostics = function(fit) list(
      num_divergent = 1L, num_max_treedepth = 0L, ebfmi = 0.8,
      unavailable = FALSE
    ),
    .package = "cpaic"
  )
  expect_warning(
    failed <- cpaic:::.cpaic_check_diagnostics(list()),
    "divergent transition"
  )
  expect_identical(failed$status, "failed")
  expect_false(failed$decision_grade)
})
