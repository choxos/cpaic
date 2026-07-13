test_that("random effects use the standard multi-arm correlation", {
  ipd <- data.frame(
    .study = rep("S1", 6),
    .trt = rep(c("A", "B", "C"), each = 2)
  )
  agd <- data.frame(
    .study = rep("S2", 2),
    .trt = c("A", "B")
  )

  re <- cpaic:::.cpaic_random_effects(ipd, agd, ".study", ".trt")
  corr <- tcrossprod(re$L_delta)

  expect_equal(re$N_delta, 3L)
  expect_equal(corr[1:2, 1:2], matrix(c(1, 0.5, 0.5, 1), 2))
  expect_equal(corr[3, 3], 1)
  expect_equal(corr[1:2, 3, drop = FALSE], matrix(0, 2, 1))
  expect_equal(re$re_idx_ipd, rep(c(0L, 1L, 2L), each = 2))
  expect_equal(re$re_idx_agd, c(0L, 3L))
})

test_that("exact survival contributions cover censoring and delayed entry", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(1, 2), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  basis <- cpaic:::.cpaic_survival_basis_eval(
    spec,
    time = rep(2, 5),
    start_time = c(0, 0, 0, 1, 0),
    entry_time = c(0, 0, 0, 0, 1)
  )
  got <- cpaic:::.cpaic_survival_loglik(
    basis, status = c(0L, 1L, 2L, 3L, 1L),
    eta = rep(log(0.2), 5), coefficients = 1
  )

  expected <- c(
    -0.4,
    log(0.2) - 0.4,
    log1p(-exp(-0.4)),
    log(exp(-0.2) - exp(-0.4)),
    log(0.2) - 0.2
  )
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("likelihood integration removes the survival mixture bias", {
  hazards <- c(0.1, 0.4)
  weights <- c(0.5, 0.5)
  follow_up <- 10
  truth <- sum(weights * (1 - exp(-hazards * follow_up)))
  mean_person_time <- sum(weights * (1 - exp(-hazards * follow_up)) / hazards)
  old <- sum(weights * hazards) * mean_person_time

  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = follow_up, baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  basis <- cpaic:::.cpaic_survival_basis_eval(
    spec, time = rep(follow_up, 2), start_time = c(0, 0),
    entry_time = c(0, 0)
  )
  node_log_lik <- cpaic:::.cpaic_survival_loglik(
    basis, status = c(2L, 2L), eta = log(hazards), coefficients = 1
  )
  exact <- sum(weights * exp(node_log_lik))

  expect_equal(truth, 0.806902460, tolerance = 1e-9)
  expect_equal(old, 1.096927061, tolerance = 1e-9)
  expect_equal(100 * (old / truth - 1), 35.942957, tolerance = 1e-6)
  expect_equal(exact, truth, tolerance = 1e-12)
})

test_that("M-spline cumulative basis differentiates to its hazard basis", {
  skip_if_not_installed("splines2")
  observed <- seq(0.1, 10, length.out = 100)
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = observed, baseline = "mspline",
    cut_points = numeric(), n_basis = 4L
  )
  eps <- 1e-5
  center <- cpaic:::.cpaic_survival_basis_eval(
    spec, time = 5, start_time = 0, entry_time = 0
  )
  left <- cpaic:::.cpaic_survival_basis_eval(
    spec, time = 5 - eps, start_time = 0, entry_time = 0
  )
  right <- cpaic:::.cpaic_survival_basis_eval(
    spec, time = 5 + eps, start_time = 0, entry_time = 0
  )

  numeric_derivative <- (right$itime - left$itime) / (2 * eps)
  expect_equal(numeric_derivative, center$time, tolerance = 1e-5)
  expect_true(all(center$time >= 0))
  expect_true(all(center$itime >= 0))
})

test_that("aggregate survival counts are rejected because they are not exact", {
  skip_if_not_installed("cmdstanr")
  ipd <- data.frame(
    .study = rep("S1", 4), .trt = rep(c("Placebo", "A"), each = 2),
    .y = c(1L, 0L, 1L, 0L), .time = c(1, 2, 1.5, 2),
    x1 = c(-1, 1, -1, 1)
  )
  agd <- data.frame(
    .study = "S2", .trt = c("Placebo", "A"),
    r = c(1L, 2L), E = c(3, 3),
    x1_mean = c(0, 0), x1_sd = c(1, 1)
  )

  expect_error(
    cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
           family = "survival"),
    "Exact survival likelihood requires `agd` event or censoring rows"
  )
})

test_that("Bayesian defaults are weakly informative and configurable", {
  f <- formals(cmlnmr)

  expect_equal(eval(f$trt_effects), c("fixed", "random"))
  expect_equal(eval(f$prior_intercept_sd), 2.5)
  expect_equal(eval(f$prior_beta_sd), 2.5)
  expect_equal(eval(f$prior_reg_sd), 1)
  expect_equal(eval(f$prior_gamma_dist), c("normal", "student_t"))
  expect_equal(eval(f$prior_gamma_scale), 1)
  expect_equal(eval(f$prior_gamma_df), 4)
  expect_equal(eval(f$prior_tau_dist), c("half-normal", "half-student-t"))
  expect_equal(eval(f$prior_tau_scale), 1)
  expect_equal(eval(f$prior_tau_df), 4)
  expect_false(eval(f$prior_predictive))
  expect_equal(eval(f$re_parameterization), c("noncentered", "centered"))
})

test_that("all Bayesian models expose pointwise log likelihood", {
  stan_dir <- system.file("stan", package = "cpaic")
  files <- list.files(stan_dir, pattern = "^cpaic_.*[.]stan$",
                      full.names = TRUE)
  expect_length(files, 5L)
  for (file in files) {
    code <- paste(readLines(file, warn = FALSE), collapse = "\n")
    expect_match(code, "vector\\[N_ipd \\+ N_agd\\] log_lik")
  }
})

test_that("LOO, WAIC, and DIC use pointwise log likelihood draws", {
  set.seed(1)
  log_lik <- matrix(rnorm(1200, -1, 0.1), ncol = 4)
  colnames(log_lik) <- paste0("log_lik[", seq_len(ncol(log_lik)), "]")
  cmdstan_fit <- list(draws = function(variables, format) log_lik)
  fit <- structure(list(fit = cmdstan_fit), class = "cpaic_mlnmr")

  expect_s3_class(loo::loo(fit), "psis_loo")
  expect_s3_class(suppressWarnings(loo::waic(fit)), "waic")
  d <- dic(fit)
  expect_s3_class(d, "cpaic_dic")
  expect_true(is.finite(d$dic))
  expect_equal(d$penalty, "pV")
})

test_that("prior sensitivity reports contrast movement from both refits", {
  C <- matrix(c(0, 1), ncol = 1,
              dimnames = list(c("Placebo", "A"), "A"))
  make_fit <- function(gamma_mean) {
    draws <- list(
      beta = matrix(rep(0.4, 100), ncol = 1,
                    dimnames = list(NULL, "beta[1]")),
      gamma = matrix(rep(gamma_mean, 100), ncol = 1,
                     dimnames = list(NULL, "gamma[1,1]"))
    )
    cmdstan_fit <- list(draws = function(variables, format) draws[[variables]])
    structure(
      list(
        fit = cmdstan_fit, C.matrix = C, effect_modifiers = "x1",
        reference = "Placebo", joint_design = diag(2),
        joint_design_ipd = diag(2),
        refit_args = list(prior_gamma_scale = 1)
      ),
      class = "cpaic_mlnmr"
    )
  }
  base <- make_fit(0.1)
  fake_refit <- function(...) {
    args <- list(...)
    make_fit(0.1 * args$prior_gamma_scale)
  }
  testthat::local_mocked_bindings(cmlnmr = fake_refit, .package = "cpaic")

  out <- prior_sensitivity(
    base, newdata = data.frame(x1 = 1), prior = "gamma",
    tighter = 0.5, looser = 2
  )

  expect_s3_class(out, "cpaic_prior_sensitivity")
  expect_equal(out$movement$estimate, 0.5)
  expect_equal(out$movement$tighter, 0.45)
  expect_equal(out$movement$looser, 0.6)
  expect_equal(out$movement$max_movement, 0.1)
})
