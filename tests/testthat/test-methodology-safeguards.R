make_safeguard_network_data <- function(family = "binomial") {
  agd <- data.frame(
    studlab = "S1", treat1 = "A", treat2 = "B",
    TE = NA_real_, seTE = NA_real_, stringsAsFactors = FALSE
  )
  ipd <- data.frame(
    .study = rep("S1", 8), .trt = rep(c("A", "B"), each = 4),
    .y = rep(c(0, 1), 4), .time = seq_len(8), .status = rep(c(0, 1), 4),
    x = seq(-1, 1, length.out = 8), stringsAsFactors = FALSE
  )
  list(agd = agd, ipd = ipd, family = family)
}

make_tied_rank_fit <- function() {
  C <- matrix(
    c(0, 0,
      1, 0,
      1, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("Placebo", "A", "A+B"), c("A", "B"))
  )
  draws <- matrix(0, nrow = 20, ncol = 2,
                  dimnames = list(NULL, c("beta[1]", "beta[2]")))
  engine <- list(draws = function(variables, format) draws)
  structure(
    list(
      fit = engine, C.matrix = C, comps = colnames(C),
      effect_modifiers = character(), reference = "Placebo",
      joint_design = diag(2), joint_design_ipd = diag(2),
      family = "binomial", sm = "OR", margins = character()
    ),
    class = c("cpaic_mlnmr", "cpaic_fit")
  )
}

test_that("target profiles reject factor and character coercion", {
  expect_error(
    cpaic:::.cpaic_target_x(data.frame(x = factor("9")), "x"),
    "numeric"
  )
  expect_error(
    cpaic:::.cpaic_target_x(data.frame(x = "9"), "x"),
    "numeric"
  )
  expect_error(
    cpaic:::.cpaic_target_x(data.frame(x = as.Date("2026-01-01")), "x"),
    "x"
  )
  expect_error(
    cpaic:::.cpaic_target_x(data.frame(x = I(list(9))), "x"),
    "x"
  )
  expect_equal(cpaic:::.cpaic_target_x(data.frame(x = 9L), "x"), c(x = 9))
})

test_that("target means respect margin support without rejecting prevalences", {
  expect_equal(
    cpaic:::.cpaic_target_x(
      data.frame(x = 0.15), "x", c(x = "bernoulli")
    ),
    c(x = 0.15)
  )
  expect_error(
    cpaic:::.cpaic_target_x(
      data.frame(x = 1.01), "x", c(x = "bernoulli")
    ),
    "[0, 1]", fixed = TRUE
  )
  expect_error(
    cpaic:::.cpaic_target_x(data.frame(x = 0), "x", c(x = "gamma")),
    "positive"
  )
  expect_error(
    cpaic:::.cpaic_target_x(data.frame(x = 1), "x", c(x = "beta")),
    "(0, 1)", fixed = TRUE
  )
})

test_that("cpaic_network rejects missing IPD values instead of row omission", {
  d <- make_safeguard_network_data()
  d$ipd$x[2] <- NA_real_
  expect_error(
    cpaic_network(d$agd, d$ipd, sm = "OR", family = "binomial",
                  ipd_covariates = "x"),
    "missing or non-finite"
  )

  d <- make_safeguard_network_data()
  d$ipd$.y[2] <- NA_real_
  expect_error(
    cpaic_network(d$agd, d$ipd, sm = "OR", family = "binomial",
                  ipd_covariates = "x"),
    "0 or 1.*no missing"
  )
})

test_that("network study identifiers cannot be missing or blank", {
  d <- make_safeguard_network_data()
  d$agd$studlab <- " "
  expect_error(
    cpaic_network(d$agd, d$ipd, sm = "OR", family = "binomial",
                  ipd_covariates = "x"),
    "missing or empty study labels"
  )

  d <- make_safeguard_network_data()
  d$ipd$.study[1] <- NA_character_
  expect_error(
    cpaic_network(d$agd, d$ipd, sm = "OR", family = "binomial",
                  ipd_covariates = "x"),
    "missing or empty labels"
  )
})

test_that("two-stage survival requires an unambiguous 0/1 event indicator", {
  d <- make_safeguard_network_data("survival")
  d$ipd$.status <- rep(c(1L, 2L), 4)
  expect_error(
    cpaic_network(
      d$agd, d$ipd, sm = "HR", family = "survival",
      ipd_time = ".time", ipd_status = ".status", ipd_covariates = "x"
    ),
    "coded 0/1"
  )
})

test_that("two-stage methods reject partial IPD from a multi-arm study", {
  ipd <- data.frame(
    .study = rep("S1", 12), .trt = rep(c("A", "B"), each = 6),
    .y = rep(c(0, 1), 6), x = seq(-1, 1, length.out = 12)
  )
  agd <- data.frame(
    studlab = rep("S1", 3),
    treat1 = c("B", "C", "C"), treat2 = c("A", "A", "B"),
    TE = c(NA, 0.3, 0.2), seTE = c(NA, 0.1, 0.1)
  )
  net <- cpaic_network(agd, ipd, sm = "OR", family = "binomial",
                       ipd_covariates = "x", reference = "A")
  expect_error(
    cstc(net, target = c(x = 0), effect_modifiers = "x"),
    "partial IPD.*multi-arm"
  )
})

test_that("ranking APIs use the same set and split exact ties", {
  fit <- make_tied_rank_fit()
  ranks <- cpaic_ranks(fit, include_screen_only = TRUE)
  probs <- rank_probs(fit, include_screen_only = TRUE)

  expect_setequal(unique(probs$element), ranks$element)
  best <- probs[probs$rank_position == 1L, c("element", "probability")]
  best <- best[match(ranks$element, best$element), ]
  expect_equal(best$probability, ranks$p_best)
  expect_equal(ranks$p_best, rep(1 / 3, 3), tolerance = 1e-12)
  prob_matrix <- xtabs(probability ~ element + rank_position, data = probs)
  expect_equal(unname(rowSums(prob_matrix)), rep(1, 3), tolerance = 1e-12)
  expect_equal(unname(colSums(prob_matrix)), rep(1, 3), tolerance = 1e-12)

  subset_probs <- rank_probs(
    fit, set = c("Placebo", "A"), include_screen_only = TRUE
  )
  expect_setequal(unique(subset_probs$element), c("Placebo", "A"))
  expect_identical(attr(ranks, "estimand"), "average_conditional_link")
  expect_identical(attr(probs, "estimand"), "average_conditional_link")
})

test_that("relative effects name link-scale uncertainty and estimand", {
  out <- relative_effects(make_tied_rank_fit())
  expect_named(out, c("treatment", "comparator", "estimate", "estimate_link",
                      "se_link", "lower", "upper", "scale", "pr_gt0",
                      "basis"))
  expect_true(all(out$scale == "natural"))
  expect_identical(attr(out, "estimand"), "average_conditional_link")
  expect_identical(attr(out, "diagnostic_status"), "unknown")
  expect_error(
    relative_effects(make_tied_rank_fit(), estimand = "marginal"),
    "does not yet implement marginal standardization", fixed = TRUE
  )
  expect_error(
    cpaic_ranks(
      make_tied_rank_fit(), include_screen_only = TRUE,
      estimand = "marginal"
    ),
    "supports estimand = \"average_conditional_link\" only", fixed = TRUE
  )
})

test_that("rank_curve preserves failed target values and their errors", {
  object <- structure(list(effect_modifiers = "x"), class = "cpaic_mlnmr")
  fake_ranks <- function(object, newdata, ...) {
    if (newdata$x == 1) stop("target-specific failure", call. = FALSE)
    data.frame(
      element = c("A", "B"), estimate = c(0, 1), p_best = c(0, 1),
      mean_rank = c(2, 1), sucra = c(0, 1)
    )
  }
  testthat::local_mocked_bindings(cpaic_ranks = fake_ranks, .package = "cpaic")

  curve <- rank_curve(object, em = "x", values = c(0, 1))
  failed <- curve[curve$status == "failed", , drop = FALSE]
  expect_equal(nrow(failed), 1L)
  expect_equal(failed$x, 1)
  expect_match(failed$error, "target-specific failure")
})

test_that("unavailable sampler diagnostics are not reported as zero", {
  testthat::local_mocked_bindings(
    .cpaic_sampler_diagnostics = function(fit) NULL,
    .cpaic_stan_variables = function(fit) "beta",
    .cpaic_fit_summary = function(fit, variables) data.frame(
      rhat = 1, ess_bulk = 200, ess_tail = 200
    ),
    .package = "cpaic"
  )
  expect_warning(
    out <- cpaic:::.cpaic_check_diagnostics(list()),
    "diagnostics.*unavailable"
  )
  expect_true(is.na(out$divergences))
  expect_true(is.na(out$max_treedepth))
})

test_that("geom_km rejects unsupported censoring and carries delayed entry", {
  skip_if_not_installed("ggplot2")
  make_fit <- function(status, entry = 0) {
    d <- data.frame(
      .study = "S1", .trt = "A", .time = c(1, 2),
      .y = status, .entry = entry
    )
    structure(
      list(
        family = "survival",
        refit_args = list(
          study = ".study", trt = ".trt", time = ".time", outcome = ".y",
          entry = ".entry", ipd = d, agd = NULL
        )
      ),
      class = c("cpaic_mlnmr", "cpaic_fit")
    )
  }

  expect_error(geom_km(make_fit(c(0L, 2L))), "right-censored")
  rows <- cpaic:::.cpaic_surv_rows(make_fit(c(0L, 1L), c(0, 0.5)))
  expect_equal(rows$entry, c(0, 0.5))
  expect_identical(unique(rows$source), "IPD")
  expect_match(unique(rows$study_arm), "^\\[IPD\\]")
  expect_length(geom_km(make_fit(c(0L, 1L), c(0, 0.5))), 2L)
})
