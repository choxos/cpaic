make_maic_bootstrap_fixture <- function(seed = 11, n = 200) {
  set.seed(seed)
  x <- rep(seq(-1, 1, length.out = n / 2), 2)
  ipd <- data.frame(
    .study = "S1",
    .trt = rep(c("P", "A"), each = n / 2),
    .y = c(stats::rnorm(n / 2), stats::rnorm(n / 2, 0.5)),
    x = x,
    stringsAsFactors = FALSE
  )
  list(
    ipd = ipd,
    info = list(trt = ".trt", outcome = ".y"),
    outcome_args = list(time = NULL, status = NULL, exposure = NULL)
  )
}

run_one_study_maic <- function(d, n_boot, min_boot_success = 0.8) {
  cpaic:::.cpaic_maic_one_study(
    d$ipd, d$info, "gaussian", "P", list(x = 0), NULL,
    "x_CENTERED", n_boot, min_boot_success, d$outcome_args, "S1"
  )
}

test_that("cMAIC retains bootstrap draws and Monte Carlo uncertainty", {
  skip_if_not_installed("maicplus")
  d <- make_maic_bootstrap_fixture()
  agd <- data.frame(
    studlab = "S1", treat1 = "A", treat2 = "P", TE = 0.5, seTE = 0.2,
    stringsAsFactors = FALSE
  )
  net <- cpaic_network(
    agd, ipd = d$ipd, sm = "MD", family = "gaussian",
    ipd_covariates = "x", inactive = "P"
  )

  fit <- cmaic(
    net, target = c(x = 0), effect_modifiers = "x",
    n_boot = 20, seed = 123
  )

  expect_true(all(c(
    "bootstrap_draws", "bootstrap_summary", "bootstrap_failures",
    "bootstrap_failure_table", "bootstrap_mcse_method",
    "bootstrap_success_rule"
  ) %in% names(fit)))
  expect_named(fit$bootstrap_draws, "S1")
  expect_equal(dim(fit$bootstrap_draws$S1), c(20L, 1L))
  expect_true(all(is.finite(fit$bootstrap_draws$S1)))
  expect_equal(fit$bootstrap_summary$n_success, 20L)
  expect_equal(fit$bootstrap_summary$success_fraction, 1)
  expect_equal(fit$bootstrap_summary$required_success_count, 20L)
  expect_equal(
    fit$bootstrap_summary$bootstrap_se_mcse_normal_approx,
    fit$bootstrap_summary$bootstrap_se / sqrt(2 * 19),
    tolerance = 1e-12
  )
  expect_match(fit$bootstrap_mcse_method, "Normal-theory Monte Carlo")
  expect_match(fit$bootstrap_success_rule, "max\\(ceiling")
  expect_equal(nrow(fit$bootstrap_failures), 0L)
  expect_equal(nrow(fit$bootstrap_failure_table), 0L)
})

test_that("every rejected cMAIC bootstrap replicate has a failure code", {
  skip_if_not_installed("maicplus")
  d <- make_maic_bootstrap_fixture()
  original_gate <- cpaic:::.cpaic_regression_problems
  gate_calls <- 0L
  mocked_gate <- function(...) {
    gate_calls <<- gate_calls + 1L
    if (gate_calls %in% c(3L, 7L)) return("GLM did not converge")
    original_gate(...)
  }

  result <- testthat::with_mocked_bindings(
    run_one_study_maic(d, n_boot = 25, min_boot_success = 0.8),
    .cpaic_regression_problems = mocked_gate,
    .package = "cpaic"
  )

  expect_equal(result$bootstrap_summary$n_success, 23L)
  expect_equal(result$bootstrap_summary$success_fraction, 23 / 25)
  expect_setequal(unique(result$bootstrap_failures$replicate), c(2L, 6L))
  expect_true(all(result$bootstrap_failures$stage == "regression_validation"))
  expect_true(all(
    result$bootstrap_failures$reason_code == "outcome_model_nonconvergence"
  ))
  expect_equal(result$bootstrap_failure_table$n_replicates, 2L)
  expect_equal(sum(!stats::complete.cases(result$bootstrap_draws)), 2L)
})

test_that("cMAIC aggregates bootstrap diagnostics across studies", {
  skip_if_not_installed("maicplus")
  d1 <- make_maic_bootstrap_fixture(n = 100)
  d2 <- make_maic_bootstrap_fixture(seed = 12, n = 100)
  d2$ipd$.study <- "S2"
  d2$ipd$.trt[d2$ipd$.trt == "A"] <- "B"
  ipd <- rbind(d1$ipd, d2$ipd)
  agd <- data.frame(
    studlab = c("S1", "S2"), treat1 = c("A", "B"), treat2 = "P",
    TE = c(0.5, 0.4), seTE = c(0.2, 0.2), stringsAsFactors = FALSE
  )
  net <- cpaic_network(
    agd, ipd = ipd, sm = "MD", family = "gaussian",
    ipd_covariates = "x", inactive = "P"
  )
  fake_one_study <- function(ipd_s, info, family, ref_arm,
                             target_mean, target_sd, em_centered_cols,
                             n_boot, min_boot_success, outcome_args,
                             study_id) {
    treatment <- setdiff(unique(as.character(ipd_s[[info$trt]])), ref_arm)
    has_failure <- identical(study_id, "S2")
    failures <- if (has_failure) {
      data.frame(
        replicate = 2L, stage = "regression_validation",
        reason_code = "outcome_model_nonconvergence",
        reason = "GLM did not converge", stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        replicate = integer(), stage = character(), reason_code = character(),
        reason = character(), stringsAsFactors = FALSE
      )
    }
    failure_table <- if (has_failure) {
      data.frame(
        stage = "regression_validation",
        reason_code = "outcome_model_nonconvergence", n_replicates = 1L,
        fraction_of_requested = 1 / n_boot, stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        stage = character(), reason_code = character(),
        n_replicates = integer(), fraction_of_requested = numeric(),
        stringsAsFactors = FALSE
      )
    }
    draws <- matrix(
      seq_len(n_boot) / 100,
      nrow = n_boot, dimnames = list(as.character(seq_len(n_boot)), treatment)
    )
    if (has_failure) draws[2, ] <- NA_real_
    n_success <- n_boot - as.integer(has_failure)
    list(
      contrasts = data.frame(
        treat1 = treatment, treat2 = ref_arm, TE = 0.5, seTE = 0.1,
        stringsAsFactors = FALSE
      ),
      ess = nrow(ipd_s),
      diagnostics = data.frame(ess = nrow(ipd_s)),
      bootstrap_draws = draws,
      bootstrap_summary = data.frame(
        treat1 = treatment, treat2 = ref_arm, n_requested = n_boot,
        n_success = n_success, success_fraction = n_success / n_boot,
        required_success_count = 20L,
        required_success_fraction = 20 / n_boot,
        bootstrap_se = 0.1,
        bootstrap_se_mcse_normal_approx = 0.01,
        stringsAsFactors = FALSE
      ),
      bootstrap_failures = failures,
      bootstrap_failure_table = failure_table,
      bootstrap_mcse_method = "Normal-theory Monte Carlo test method",
      bootstrap_success_rule = "test success rule"
    )
  }

  fit <- testthat::with_mocked_bindings(
    suppressWarnings(cmaic(
      net, target = c(x = 0), effect_modifiers = "x",
      n_boot = 25, seed = 9
    )),
    .cpaic_maic_one_study = fake_one_study,
    .package = "cpaic"
  )

  expect_setequal(names(fit$bootstrap_draws), c("S1", "S2"))
  expect_equal(vapply(fit$bootstrap_draws, nrow, integer(1)), c(S1 = 25L, S2 = 25L))
  expect_setequal(fit$bootstrap_summary$study, c("S1", "S2"))
  expect_identical(fit$bootstrap_failures$study, "S2")
  expect_identical(fit$bootstrap_failure_table$study, "S2")
  expect_equal(fit$bootstrap_failure_table$n_replicates, 1L)
})

test_that("cMAIC bootstrap floor retains failed-fit diagnostics", {
  skip_if_not_installed("maicplus")
  d <- make_maic_bootstrap_fixture()
  original_gate <- cpaic:::.cpaic_regression_problems
  gate_calls <- 0L
  mocked_gate <- function(...) {
    gate_calls <<- gate_calls + 1L
    if (gate_calls == 2L) return("GLM did not converge")
    original_gate(...)
  }

  condition <- tryCatch(
    testthat::with_mocked_bindings(
      run_one_study_maic(d, n_boot = 10, min_boot_success = 0.5),
      .cpaic_regression_problems = mocked_gate,
      .package = "cpaic"
    ),
    cpaic_bootstrap_error = identity
  )

  expect_s3_class(condition, "cpaic_bootstrap_error")
  expect_match(conditionMessage(condition), "required 10")
  expect_equal(condition$bootstrap_summary$n_success, 9L)
  expect_equal(condition$bootstrap_summary$required_success_count, 10L)
  expect_match(condition$bootstrap_success_rule, "all.*must succeed")
  expect_equal(unique(condition$bootstrap_failures$replicate), 1L)
  expect_equal(
    condition$bootstrap_failures$reason_code,
    "outcome_model_nonconvergence"
  )
})
