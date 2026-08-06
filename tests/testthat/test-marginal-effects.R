make_marginal_draw_matrix <- function(x, names) {
  out <- if (is.matrix(x)) x else matrix(as.numeric(x), ncol = 1L)
  stopifnot(ncol(out) == length(names))
  colnames(out) <- names
  out
}

make_marginal_fit <- function(
    family = "binomial", beta = rep(log(2), 8L),
    gamma = rep(log(1.5), length(beta)),
    breg = rep(0.8, length(beta)),
    mu = matrix(rep(-1, length(beta)), ncol = 1L),
    study_levels = "S1", joint_design = diag(2),
    joint_design_ipd = joint_design, coefficients = NULL,
    margins = c(x = "normal"), survival_spec = NULL) {
  beta <- make_marginal_draw_matrix(beta, "beta[1]")
  gamma <- make_marginal_draw_matrix(gamma, "gamma[1,1]")
  breg <- make_marginal_draw_matrix(breg, "breg[1]")
  mu <- make_marginal_draw_matrix(
    mu, paste0("mu[", seq_along(study_levels), "]")
  )
  blocks <- list(beta = beta, gamma = gamma, breg = breg, mu = mu)
  if (!is.null(coefficients)) {
    blocks$coefficients <- make_marginal_draw_matrix(
      coefficients,
      unlist(lapply(seq_along(study_levels), function(s) {
        paste0("coefficients[", s, ",",
               seq_len(ncol(coefficients) / length(study_levels)), "]")
      }))
    )
  }
  stopifnot(length(unique(vapply(blocks, nrow, integer(1)))) == 1L)
  engine <- list(draws = function(variables, format = "draws_matrix") {
    if (length(variables) == 1L) {
      out <- blocks[[variables]]
      if (is.null(out)) stop("unknown draw block: ", variables)
      return(out)
    }
    do.call(cbind, blocks[variables])
  })
  C <- matrix(c(0, 1), ncol = 1L,
              dimnames = list(c("Placebo", "A"), "A"))
  structure(
    list(
      fit = engine, C.matrix = C, comps = "A",
      family = family, effect_modifiers = "x", margins = margins,
      cor = matrix(1, 1, 1, dimnames = list("x", "x")),
      reference = "Placebo",
      sm = switch(family, binomial = "OR", gaussian = "MD",
                  poisson = "IRR", survival = "HR"),
      joint_design = joint_design, joint_design_ipd = joint_design_ipd,
      study_levels = study_levels, survival_spec = survival_spec,
      survival_study_support = if (is.null(survival_spec)) NULL else
        stats::setNames(rep(survival_spec$max_time, length(study_levels)),
                        study_levels),
      trt_effects = "fixed", diagnostics = list(status = "passed"),
      refit_args = NULL, observed = NULL
    ),
    class = c("cpaic_mlnmr", "cpaic_fit")
  )
}

binomial_marginal_truth <- function(mu, beta, breg, gamma, x, weights,
                                    measure = "odds_ratio") {
  w <- weights / sum(weights)
  p0 <- sum(w * stats::plogis(mu + breg * x))
  p1 <- sum(w * stats::plogis(mu + breg * x + beta + gamma * x))
  switch(
    measure,
    odds_ratio = (p1 / (1 - p1)) / (p0 / (1 - p0)),
    risk_ratio = p1 / p0,
    risk_difference = p1 - p0
  )
}

poisson_marginal_truth <- function(beta, breg, gamma, x, weights,
                                   measure = "rate_ratio", mu = 0) {
  w <- weights / sum(weights)
  r0 <- sum(w * exp(mu + breg * x))
  r1 <- sum(w * exp(mu + breg * x + beta + gamma * x))
  switch(measure, rate_ratio = r1 / r0, rate_difference = r1 - r0)
}

test_that("gaussian marginal mean differences have the identity-link truth", {
  fit <- make_marginal_fit(
    family = "gaussian", beta = rep(0.7, 8), gamma = rep(0.4, 8),
    breg = rep(4, 8), mu = matrix(rep(-5, 8), ncol = 1)
  )
  target <- data.frame(x = c(-1, 2))
  weights <- c(0.4, 0.6)
  truth <- 0.7 + 0.4 * weighted.mean(target$x, weights)

  out <- marginal_effects(fit, target, weights = weights)
  expect_s3_class(out, "cpaic_effects")
  expect_equal(out$estimate, truth, tolerance = 1e-12)
  expect_equal(out$estimate_contrast, truth, tolerance = 1e-12)
  expect_identical(out$contrast_scale, "natural difference")
  expect_identical(attr(out, "estimand"), "marginal")
  expect_identical(attr(out, "measure"), "mean_difference")
})

test_that("binomial marginal measures match response-scale standardization", {
  fit <- make_marginal_fit()
  target <- data.frame(x = c(0, 1))
  weights <- c(0.25, 0.75)
  truth <- vapply(
    c("odds_ratio", "risk_ratio", "risk_difference"),
    function(m) binomial_marginal_truth(
      -1, log(2), 0.8, log(1.5), target$x, weights, m
    ), numeric(1)
  )

  for (m in names(truth)) {
    out <- marginal_effects(
      fit, target, weights = weights, measure = m, baseline_study = "S1"
    )
    expect_equal(out$estimate, unname(truth[[m]]), tolerance = 1e-12,
                 info = m)
    expected_link <- if (m == "risk_difference") truth[[m]] else log(truth[[m]])
    expect_equal(out$estimate_contrast, unname(expected_link),
                 tolerance = 1e-12, info = m)
  }
  expect_equal(
    marginal_effects(fit, target, weights = weights,
                     baseline_study = "S1")$estimate,
    unname(truth[["odds_ratio"]]), tolerance = 1e-12
  )
})

test_that("binomial marginal effects use the selected baseline study", {
  mu <- cbind(rep(-3, 8), rep(1, 8))
  fit <- make_marginal_fit(mu = mu, study_levels = c("S-low", "S-high"))
  target <- data.frame(x = c(0, 1))
  weights <- c(0.25, 0.75)
  low <- marginal_effects(fit, target, weights = weights,
                          baseline_study = "S-low")
  high <- marginal_effects(fit, target, weights = weights,
                           baseline_study = "S-high")

  expect_equal(
    low$estimate,
    binomial_marginal_truth(-3, log(2), 0.8, log(1.5), target$x, weights),
    tolerance = 1e-12
  )
  expect_equal(
    high$estimate,
    binomial_marginal_truth(1, log(2), 0.8, log(1.5), target$x, weights),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(low$estimate, high$estimate)))
})

test_that("poisson marginal rate ratios include prognostic effects", {
  fit <- make_marginal_fit(
    family = "poisson", beta = rep(log(1.8), 8), gamma = rep(0.2, 8),
    breg = rep(0.3, 8), mu = matrix(rep(6, 8), ncol = 1)
  )
  target <- data.frame(x = c(-1, 2))
  weights <- c(0.4, 0.6)
  truth <- poisson_marginal_truth(log(1.8), 0.3, 0.2, target$x, weights)
  out <- marginal_effects(fit, target, weights = weights)

  expect_equal(out$estimate, truth, tolerance = 1e-12)
  expect_equal(out$estimate_contrast, log(truth), tolerance = 1e-12)
  expect_identical(out$contrast_scale, "log ratio")
  expect_identical(attr(out, "measure"), "rate_ratio")
})

test_that("poisson rate differences use an explicit donor baseline", {
  fit <- make_marginal_fit(
    family = "poisson", beta = rep(log(1.8), 8), gamma = rep(0.2, 8),
    breg = rep(0.3, 8), mu = matrix(rep(log(0.2), 8), ncol = 1)
  )
  target <- data.frame(x = c(-1, 2))
  weights <- c(0.4, 0.6)
  truth <- poisson_marginal_truth(
    log(1.8), 0.3, 0.2, target$x, weights,
    measure = "rate_difference", mu = log(0.2)
  )
  out <- marginal_effects(
    fit, target, weights = weights, measure = "rate_difference",
    baseline_study = "S1"
  )
  expect_equal(out$estimate, truth, tolerance = 1e-12)
  expect_equal(out$estimate_contrast, truth, tolerance = 1e-12)
  expect_identical(out$contrast_scale, "natural difference")
  expect_error(
    marginal_effects(fit, target, weights = weights,
                     measure = "rate_difference"),
    "baseline_study"
  )
})

test_that("marginal summaries standardize each posterior draw before averaging", {
  mu <- matrix(c(-2, 1), ncol = 1)
  beta <- log(c(1.5, 3))
  fit <- make_marginal_fit(
    beta = beta, gamma = c(0.2, 0.2), breg = c(0.4, 0.4), mu = mu
  )
  target <- data.frame(x = c(-1, 1))
  weights <- c(0.3, 0.7)
  draw_truth <- vapply(seq_along(beta), function(i) {
    binomial_marginal_truth(
      mu[i, 1], beta[i], 0.4, 0.2, target$x, weights
    )
  }, numeric(1))
  plug_in <- binomial_marginal_truth(
    mean(mu), mean(beta), 0.4, 0.2, target$x, weights
  )
  out <- marginal_effects(fit, target, weights = weights,
                          baseline_study = "S1")

  expect_equal(out$estimate, mean(draw_truth), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(out$estimate, plug_in)))
})

test_that("marginal effects depend on the target distribution, not only its mean", {
  point <- data.frame(x = 0)
  spread <- data.frame(x = c(-1, 1))
  w <- c(0.5, 0.5)

  bf <- make_marginal_fit(
    mu = matrix(rep(-0.5, 8), ncol = 1), beta = rep(log(2), 8),
    gamma = rep(log(1.5), 8), breg = rep(0.8, 8)
  )
  b0 <- marginal_effects(bf, point, baseline_study = "S1")
  bs <- marginal_effects(bf, spread, weights = w, baseline_study = "S1")
  expect_equal(b0$estimate, 2, tolerance = 1e-12)
  expect_equal(bs$estimate, 1.762004204794788, tolerance = 1e-12)

  pf <- make_marginal_fit(
    family = "poisson", beta = rep(log(1.5), 8),
    gamma = rep(log(1.4), 8), breg = rep(0.6, 8)
  )
  p0 <- marginal_effects(pf, point)
  ps <- marginal_effects(pf, spread, weights = w)
  expect_equal(p0$estimate, 1.5, tolerance = 1e-12)
  expect_equal(ps$estimate, 1.861911205884704, tolerance = 1e-12)

  gf <- make_marginal_fit(
    family = "gaussian", beta = rep(0.7, 8), gamma = rep(0.4, 8)
  )
  expect_equal(marginal_effects(gf, point)$estimate,
               marginal_effects(gf, spread, weights = w)$estimate,
               tolerance = 1e-12)
})

test_that("empirical target integration is permutation and weight invariant", {
  fit <- make_marginal_fit(family = "poisson")
  target <- data.frame(x = c(-1, 0.5, 2))
  weights <- c(1, 2, 4)
  a <- marginal_effects(fit, target, weights = weights)
  b <- marginal_effects(fit, target[3:1, , drop = FALSE],
                        weights = weights[3:1])
  c <- marginal_effects(fit, target, weights = 100 * weights)
  columns <- function(x) stats::setNames(
    lapply(names(x), function(nm) x[[nm]]), names(x)
  )

  expect_equal(columns(a), columns(b))
  expect_equal(columns(a), columns(c))
  expect_equal(attr(a, "target_mean"), attr(b, "target_mean"))
  expect_equal(
    columns(marginal_effects(fit, target)),
    columns(marginal_effects(
      fit, target, weights = rep(1, nrow(target))
    ))
  )
})

test_that("summary targets reproduce their explicit integration nodes", {
  skip_if_not_installed("randtoolbox")
  fit <- make_marginal_fit(family = "poisson")
  spec <- list(
    means = c(x = 0.2), sds = c(x = 1.1),
    margins = c(x = "normal"), n_int = 64L
  )
  summary_result <- marginal_effects(fit, spec)
  resolved <- attr(summary_result, "target_distribution")
  node_result <- marginal_effects(
    fit, resolved$nodes, weights = resolved$weights
  )
  columns <- function(x) stats::setNames(
    lapply(names(x), function(nm) x[[nm]]), names(x)
  )
  expect_equal(columns(summary_result), columns(node_result))
  expect_equal(attr(summary_result, "target_mean"),
               attr(node_result, "target_mean"))
  expect_identical(attr(summary_result, "target_distribution")$type,
                   "summary")
  expect_identical(attr(node_result, "target_distribution")$type,
                   "empirical")
  expect_equal(attr(summary_result, "target_mean"), spec$means,
               tolerance = 1e-12)
})

test_that("summary targets preserve Gaussian means and rare Bernoulli strata", {
  normal_fit <- make_marginal_fit(
    family = "gaussian", beta = rep(0.7, 8), gamma = rep(0.4, 8)
  )
  normal <- marginal_effects(
    normal_fit,
    list(means = c(x = 0.2), sds = c(x = 1.1), n_int = 2L)
  )
  expect_equal(normal$estimate, 0.7 + 0.4 * 0.2, tolerance = 1e-12)
  expect_equal(attr(normal, "target_mean"), c(x = 0.2), tolerance = 1e-12)

  rare_fit <- make_marginal_fit(
    family = "poisson", margins = c(x = "bernoulli")
  )
  rare <- marginal_effects(
    rare_fit,
    list(means = c(x = 0.001), sds = c(x = 0),
         margins = c(x = "bernoulli"), n_int = 1L)
  )
  resolved <- attr(rare, "target_distribution")
  expect_identical(sort(resolved$nodes$x), c(0, 1))
  expect_equal(sum(resolved$weights * resolved$nodes$x), 0.001,
               tolerance = 1e-15)
})

test_that("marginal target and argument validation fail at the front door", {
  bf <- make_marginal_fit()
  pf <- make_marginal_fit(family = "poisson")
  expect_error(marginal_effects(list(), data.frame(x = 0)), "cmlnmr")
  expect_error(marginal_effects(pf, NULL), "target")
  expect_error(marginal_effects(pf, data.frame(y = 0)), "missing.*x")
  expect_error(marginal_effects(pf, data.frame(x = 0, y = 1)), "exactly")
  expect_error(marginal_effects(pf, data.frame(x = factor(0))), "numeric")
  expect_error(marginal_effects(pf, data.frame(x = "0")), "numeric")
  expect_error(marginal_effects(pf, data.frame(x = NA_real_)), "finite")
  expect_error(marginal_effects(pf, data.frame(x = Inf)), "finite")
  expect_error(marginal_effects(pf, data.frame(x = numeric())), "row")

  target <- data.frame(x = c(0, 1))
  expect_error(marginal_effects(pf, target, weights = 1), "weights")
  expect_error(marginal_effects(pf, target, weights = c(1, -1)), "weights")
  expect_error(marginal_effects(pf, target, weights = c(0, 0)), "weights")
  expect_error(marginal_effects(pf, target, weights = c(1, NA)), "weights")
  expect_error(
    marginal_effects(pf, list(means = c(x = 0), sds = c(x = 1)),
                     weights = 1),
    "weights"
  )
  expect_error(marginal_effects(pf, target, measure = "odds_ratio"),
               "measure")
  expect_error(marginal_effects(bf, target, baseline_study = "S1",
                                measure = "rate_ratio"), "measure")
  expect_error(marginal_effects(bf, target), "baseline_study")
  expect_error(marginal_effects(bf, target, baseline_study = "unknown"),
               "baseline_study")
  expect_error(marginal_effects(pf, target, reference = "unknown"),
               "reference")
  expect_error(marginal_effects(pf, target, all_contrasts = NA),
               "TRUE or FALSE")
  expect_error(marginal_effects(pf, target, level = 1), "\\(0, 1\\)")
  expect_error(marginal_effects(pf, target, times = 1), "only for survival")
  expect_error(marginal_effects(pf, target, baseline_study = "S1"),
               "not used")
  expect_error(marginal_effects(pf, target, typo = TRUE), "Unused argument")

  bad_cor <- matrix(1, nrow = 1, dimnames = list("wrong", "wrong"))
  expect_error(
    marginal_effects(
      pf, list(means = c(x = 0), sds = c(x = 1), cor = bad_cor)
    ),
    "dimnames"
  )
})

test_that("marginal prediction extracts posterior blocks jointly", {
  fit <- make_marginal_fit(family = "poisson")
  original_draws <- fit$fit$draws
  calls <- list()
  fit$fit$draws <- function(variables, format = "draws_matrix") {
    calls[[length(calls) + 1L]] <<- variables
    original_draws(variables, format)
  }
  expect_s3_class(marginal_effects(fit, data.frame(x = 0)), "cpaic_effects")
  expect_length(calls, 1L)
  expect_setequal(calls[[1L]], c("beta", "breg", "gamma", "mu"))

  incomplete <- make_marginal_fit(family = "poisson")
  original_draws <- incomplete$fit$draws
  incomplete$fit$draws <- function(variables, format = "draws_matrix") {
    out <- original_draws(variables, format)
    out[, !grepl("^gamma\\[", colnames(out)), drop = FALSE]
  }
  expect_error(
    marginal_effects(incomplete, data.frame(x = 0)),
    "component x modifier"
  )
})

test_that("empirical target weights normalize without overflow", {
  fit <- make_marginal_fit(family = "poisson")
  target <- data.frame(x = c(-1, 1))
  ordinary <- marginal_effects(fit, target, weights = c(1, 1))
  huge <- marginal_effects(
    fit, target, weights = rep(.Machine$double.xmax, 2L)
  )
  expect_equal(huge$estimate, ordinary$estimate)
  expect_equal(huge$estimate_contrast, ordinary$estimate_contrast)
  expect_equal(attr(huge, "target_distribution")$weights, c(0.5, 0.5))
})

test_that("correlated summary targets fail closed on joint distortion", {
  fit <- make_marginal_fit(family = "poisson")
  fit$effect_modifiers <- c("b", "x")
  fit$margins <- c(b = "bernoulli", x = "normal")
  target_cor <- matrix(
    c(1, 0.9, 0.9, 1), nrow = 2L,
    dimnames = list(c("b", "x"), c("b", "x"))
  )
  fit$cor <- target_cor
  target <- list(
    means = c(b = 0.5, x = 0), sds = c(b = 0.5, x = 1),
    margins = fit$margins, cor = target_cor, n_int = 2L
  )
  expect_error(marginal_effects(fit, target), "joint dependence")
})

test_that("empirical targets enforce fitted marginal support", {
  bern <- make_marginal_fit(family = "poisson", margins = c(x = "bernoulli"))
  expect_error(marginal_effects(bern, data.frame(x = 0.25)), "0 or 1")
  expect_s3_class(
    marginal_effects(bern, data.frame(x = c(0, 1))), "cpaic_effects"
  )
  gamma <- make_marginal_fit(family = "poisson", margins = c(x = "gamma"))
  expect_error(marginal_effects(gamma, data.frame(x = c(0, 1))), "positive")
  beta <- make_marginal_fit(family = "poisson", margins = c(x = "beta"))
  expect_error(marginal_effects(beta, data.frame(x = c(0, 0.5))),
               "\\(0, 1\\)")
})

test_that("marginal estimability uses target support for nonlinear links", {
  beta_only <- matrix(c(1, 0), nrow = 1L)
  target <- data.frame(x = c(-1, 1))
  weights <- c(0.5, 0.5)

  for (family in c("binomial", "poisson")) {
    fit <- make_marginal_fit(
      family = family, joint_design = beta_only,
      joint_design_ipd = beta_only
    )
    args <- list(object = fit, target = target, weights = weights)
    if (family == "binomial") args$baseline_study <- "S1"
    out <- do.call(marginal_effects, args)
    expect_true(is.na(out$estimate), info = family)
    expect_identical(out$basis, "first-order screen failed", info = family)
  }

  gaussian <- make_marginal_fit(
    family = "gaussian", joint_design = beta_only,
    joint_design_ipd = beta_only
  )
  gout <- marginal_effects(gaussian, target, weights = weights)
  expect_true(is.finite(gout$estimate))
  expect_identical(gout$basis, "exact")
})

test_that("nonlinear marginal identification is labeled as a screen", {
  target <- data.frame(x = 1)
  aggregate_gamma <- matrix(c(1, 0), nrow = 1L)
  fit <- make_marginal_fit(
    family = "poisson", joint_design = diag(2),
    joint_design_ipd = aggregate_gamma
  )
  out <- marginal_effects(fit, target)
  expect_true(is.finite(out$estimate))
  expect_identical(out$basis, "first-order screen")
})

test_that("marginal effects retain study metadata after redaction", {
  fit <- make_marginal_fit()
  redacted <- redact_fit(fit)
  expect_null(redacted$refit_args)
  expect_identical(redacted$study_levels, "S1")
  expect_equal(
    marginal_effects(redacted, data.frame(x = 0), baseline_study = "S1"),
    marginal_effects(fit, data.frame(x = 0), baseline_study = "S1")
  )

  legacy <- fit
  legacy$study_levels <- NULL
  expect_error(
    marginal_effects(legacy, data.frame(x = 0), baseline_study = "S1"),
    "study.*metadata"
  )
})

test_that("marginal calculations are stable at extreme predictors", {
  bin <- make_marginal_fit(
    beta = rep(log(2), 8), gamma = rep(0, 8), breg = rep(0, 8),
    mu = matrix(rep(50, 8), ncol = 1)
  )
  bout <- marginal_effects(bin, data.frame(x = 0), baseline_study = "S1")
  expect_true(is.finite(bout$estimate))
  expect_equal(bout$estimate, 2, tolerance = 1e-10)

  pois <- make_marginal_fit(
    family = "poisson", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(1000, 8)
  )
  pout <- marginal_effects(pois, data.frame(x = 1))
  expect_true(is.finite(pout$estimate))
  expect_equal(pout$estimate, 2, tolerance = 1e-10)
})

test_that("relative_effects marginal dispatch matches the dedicated API", {
  fit <- make_marginal_fit()
  target <- data.frame(x = c(0, 1))
  weights <- c(1, 3)
  direct <- marginal_effects(
    fit, target, weights = weights, measure = "risk_ratio",
    baseline_study = "S1"
  )
  wrapped <- relative_effects(
    fit, estimand = "marginal", target = target, weights = weights,
    measure = "risk_ratio", baseline_study = "S1"
  )
  expect_equal(wrapped, direct)
  expect_identical(attr(wrapped, "estimand"), "marginal")
  expect_identical(attr(wrapped, "measure"), "risk_ratio")
})

test_that("backtransf controls the reported marginal scale", {
  fit <- make_marginal_fit()
  target <- data.frame(x = c(0, 1))
  natural <- marginal_effects(
    fit, target, baseline_study = "S1", measure = "odds_ratio"
  )
  link <- marginal_effects(
    fit, target, baseline_study = "S1", measure = "odds_ratio",
    backtransf = FALSE
  )
  expect_equal(link$estimate, natural$estimate_contrast, tolerance = 1e-12)
  expect_equal(link$estimate, log(natural$estimate), tolerance = 1e-12)
  expect_true(all(link$scale == "link"))
})

make_piecewise_marginal_survival_fit <- function() {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(1, 10), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  spec$max_time <- 10
  make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(log(0.1), 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = spec
  )
}

test_that("piecewise survival and risk contrasts have closed-form truths", {
  fit <- make_piecewise_marginal_survival_fit()
  times <- c(2, 5)
  target <- data.frame(x = 0)
  s0 <- exp(-0.1 * times)
  s1 <- exp(-0.2 * times)
  truths <- list(
    survival_difference = s1 - s0,
    survival_ratio = s1 / s0,
    risk_difference = (1 - s1) - (1 - s0),
    risk_ratio = (1 - s1) / (1 - s0),
    time_specific_hazard_ratio = rep(2, length(times))
  )

  for (measure in names(truths)) {
    out <- marginal_effects(
      fit, target, measure = measure, baseline_study = "S1", times = times
    )
    expect_equal(out$time, times, tolerance = 0, info = measure)
    expect_equal(out$estimate, truths[[measure]], tolerance = 1e-12,
                 info = measure)
  }
})

test_that("piecewise RMST contrasts have closed-form truths", {
  fit <- make_piecewise_marginal_survival_fit()
  times <- c(2, 5)
  target <- data.frame(x = 0)
  rmst <- function(rate, horizon) -expm1(-rate * horizon) / rate
  r0 <- rmst(0.1, times)
  r1 <- rmst(0.2, times)

  difference <- marginal_effects(
    fit, target, measure = "rmst_difference", baseline_study = "S1",
    times = times
  )
  ratio <- marginal_effects(
    fit, target, measure = "rmst_ratio", baseline_study = "S1",
    times = times
  )
  expect_equal(difference$time, times)
  expect_equal(difference$estimate, r1 - r0, tolerance = 1e-10)
  expect_equal(ratio$estimate, r1 / r0, tolerance = 1e-10)
})

test_that("survival ratios retain finite log contrasts beyond response range", {
  high <- make_piecewise_marginal_survival_fit()
  high$fit <- make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(log(1000), 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = high$survival_spec
  )$fit
  high$survival_study_support <- c(S1 = 10)
  log_survival_ratio <- marginal_effects(
    high, data.frame(x = 0), measure = "survival_ratio",
    baseline_study = "S1", times = 1, backtransf = FALSE
  )
  expect_true(is.finite(log_survival_ratio$estimate))
  expect_equal(log_survival_ratio$estimate, -1000, tolerance = 1e-10)

  low <- make_piecewise_marginal_survival_fit()
  low$fit <- make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(-1000, 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = low$survival_spec
  )$fit
  low$survival_study_support <- c(S1 = 10)
  log_risk_ratio <- marginal_effects(
    low, data.frame(x = 0), measure = "risk_ratio",
    baseline_study = "S1", times = 1, backtransf = FALSE
  )
  expect_true(is.finite(log_risk_ratio$estimate))
  expect_equal(log_risk_ratio$estimate, log(2), tolerance = 1e-10)
})

test_that("risk differences retain the complementary survival information", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(0.1, 2), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  spec$max_time <- 2
  fit <- make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(4, 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = spec
  )
  truth <- exp(-exp(4)) - exp(-2 * exp(4))
  risk <- marginal_effects(
    fit, data.frame(x = 0), measure = "risk_difference",
    baseline_study = "S1", times = 1
  )
  survival <- marginal_effects(
    fit, data.frame(x = 0), measure = "survival_difference",
    baseline_study = "S1", times = 1
  )
  expect_equal(risk$estimate, truth, tolerance = 1e-12)
  expect_equal(survival$estimate, -truth, tolerance = 1e-12)
})

test_that("survival marginal effects require a donor baseline and valid times", {
  fit <- make_piecewise_marginal_survival_fit()
  target <- data.frame(x = 0)
  expect_error(marginal_effects(fit, target), "explicit.*measure")
  expect_error(
    marginal_effects(fit, target, measure = "survival_difference", times = 5),
    "baseline_study"
  )
  expect_error(
    marginal_effects(fit, target, measure = "survival_difference",
                     baseline_study = "unknown", times = 5),
    "baseline_study"
  )
  expect_error(
    marginal_effects(fit, target, measure = "survival_difference",
                     baseline_study = "S1"),
    "times"
  )
  for (bad in list(0, -1, NA_real_, Inf, 10.01)) {
    expect_error(
      marginal_effects(fit, target, measure = "survival_difference",
                       baseline_study = "S1", times = bad),
      "times|support|10", info = paste(bad, collapse = ",")
    )
  }
  expect_s3_class(
    marginal_effects(fit, target, measure = "survival_difference",
                     baseline_study = "S1", times = 10),
    "cpaic_effects"
  )
})

test_that("survival prediction respects the selected donor's follow-up", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(1, 10), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  fit <- make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8),
    mu = cbind(rep(log(0.1), 8), rep(log(0.1), 8)),
    study_levels = c("S-short", "S-long"),
    coefficients = matrix(1, nrow = 8, ncol = 2),
    survival_spec = spec
  )
  fit$survival_study_support <- c(`S-short` = 5, `S-long` = 10)

  expect_error(
    marginal_effects(
      fit, data.frame(x = 0), measure = "survival_difference",
      baseline_study = "S-short", times = 6
    ),
    "S-short.*5|support"
  )
  expect_s3_class(
    marginal_effects(
      fit, data.frame(x = 0), measure = "survival_difference",
      baseline_study = "S-long", times = 6
    ),
    "cpaic_effects"
  )

  fit$survival_study_support <- NULL
  expect_error(
    marginal_effects(
      fit, data.frame(x = 0), measure = "survival_difference",
      baseline_study = "S-long", times = 5
    ),
    "follow-up support"
  )
})

test_that("survival rejects a scalar marginal hazard ratio", {
  fit <- make_piecewise_marginal_survival_fit()
  expect_error(
    marginal_effects(
      fit, data.frame(x = 0), measure = "hazard_ratio",
      baseline_study = "S1", times = 5
    ),
    "time_specific_hazard_ratio|measure"
  )
})

test_that("survival wrapper parity uses retained basis metadata", {
  fit <- make_piecewise_marginal_survival_fit()
  fit <- redact_fit(fit)
  expect_null(fit$refit_args)
  expect_s3_class(fit$survival_spec, "cpaic_survival_basis")
  direct <- marginal_effects(
    fit, data.frame(x = 0), measure = "risk_ratio",
    baseline_study = "S1", times = c(2, 5)
  )
  wrapped <- relative_effects(
    fit, estimand = "marginal", target = data.frame(x = 0),
    measure = "risk_ratio", baseline_study = "S1", times = c(2, 5)
  )
  expect_equal(wrapped, direct)
})

test_that("survival posterior blocks are extracted in one aligned draw set", {
  fit <- make_piecewise_marginal_survival_fit()
  direct <- marginal_effects(
    fit, data.frame(x = c(-1, 1)), weights = c(0.4, 0.6),
    measure = "survival_difference", baseline_study = "S1", times = 4
  )
  original_draws <- fit$fit$draws
  calls <- list()
  fit$fit$draws <- function(variables, format = "draws_matrix") {
    calls[[length(calls) + 1L]] <<- variables
    out <- original_draws(variables, format)
    if (identical(variables, "coefficients")) {
      out <- out[nrow(out):1L, , drop = FALSE]
    }
    out
  }
  aligned <- marginal_effects(
    fit, data.frame(x = c(-1, 1)), weights = c(0.4, 0.6),
    measure = "survival_difference", baseline_study = "S1", times = 4
  )
  expect_equal(aligned$estimate, direct$estimate)
  expect_length(calls, 1L)
  expect_setequal(
    calls[[1L]], c("beta", "breg", "mu", "gamma", "coefficients")
  )
})

test_that("marginal hazard ratios vary over time under prognostic heterogeneity", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(0.1, 10), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  spec$max_time <- 10
  fit <- make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(log(4), 8), mu = matrix(rep(log(0.1), 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = spec
  )
  target <- data.frame(x = c(0, 1))
  weights <- c(0.5, 0.5)
  times <- c(0.5, 5)

  # The conditional treatment HR is exactly two in both covariate strata.
  conditional_hr <- exp(log(2) + 0 * target$x)
  expect_equal(conditional_hr, c(2, 2))

  # Independently standardize the stratum-specific hazards over each risk set.
  marginal_hazard <- function(time, treatment_hr) {
    individual_hazard <- 0.1 * treatment_hr * exp(log(4) * target$x)
    individual_survival <- exp(-individual_hazard * time)
    sum(weights * individual_hazard * individual_survival) /
      sum(weights * individual_survival)
  }
  truth <- vapply(times, function(time) {
    marginal_hazard(time, 2) / marginal_hazard(time, 1)
  }, numeric(1))

  out <- marginal_effects(
    fit, target, weights = weights,
    measure = "time_specific_hazard_ratio",
    baseline_study = "S1", times = times
  )
  expect_equal(out$time, times)
  expect_equal(out$estimate, truth, tolerance = 1e-12)
  expect_gt(abs(diff(out$estimate)), 0.1)
  expect_false(any(abs(out$estimate - conditional_hr[1]) < 1e-8))
})

test_that("marginal hazard ratios cancel enormous common cumulative hazards", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(0.1, 2), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  spec$max_time <- 2
  for (common_eta in c(30, 35, 36, 37, 40, 50, 710)) {
    fit <- make_marginal_fit(
      family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
      breg = rep(0, 8),
      mu = matrix(rep(common_eta, 8), ncol = 1),
      coefficients = matrix(1, nrow = 8, ncol = 1),
      survival_spec = spec
    )
    out <- marginal_effects(
      fit, data.frame(x = 0),
      measure = "time_specific_hazard_ratio",
      baseline_study = "S1", times = 1
    )
    expect_equal(out$estimate_contrast, log(2), tolerance = 1e-12)
    expect_equal(out$estimate, 2, tolerance = 1e-12)
  }
})

test_that("survival ratios preserve exact equality after survival underflow", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(0.1, 2), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  spec$max_time <- 2
  fit <- make_marginal_fit(
    family = "survival", beta = rep(0, 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(710, 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = spec
  )
  out <- marginal_effects(
    fit, data.frame(x = 0), measure = "survival_ratio",
    baseline_study = "S1", times = 1
  )
  expect_equal(out$estimate_contrast, 0, tolerance = 0)
  expect_equal(out$estimate, 1, tolerance = 0)
})

test_that("piecewise RMST ratios retain subnormal information", {
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = c(0.1, 2), baseline = "piecewise",
    cut_points = numeric(), n_basis = 1L
  )
  spec$max_time <- 2
  fit <- make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(710, 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 1),
    survival_spec = spec
  )
  out <- marginal_effects(
    fit, data.frame(x = 0), measure = "rmst_ratio",
    baseline_study = "S1", times = 1, backtransf = FALSE
  )
  expect_equal(out$estimate, -log(2), tolerance = 1e-12)
  expect_equal(out$estimate_contrast, -log(2), tolerance = 1e-12)
})

make_mspline_marginal_survival_fit <- function() {
  skip_if_not_installed("splines2")
  spec <- cpaic:::.cpaic_survival_basis_spec(
    observed_times = seq(0.1, 2, length.out = 30),
    baseline = "mspline", n_basis = 4L
  )
  spec$max_time <- 2
  make_marginal_fit(
    family = "survival", beta = rep(log(2), 8), gamma = rep(0, 8),
    breg = rep(0, 8), mu = matrix(rep(log(0.1), 8), ncol = 1),
    coefficients = matrix(1, nrow = 8, ncol = 4),
    survival_spec = spec
  )
}

test_that("M-spline marginal RMST differences are rejected", {
  fit <- make_mspline_marginal_survival_fit()
  expect_error(
    marginal_effects(
      fit, data.frame(x = 0), measure = "rmst_difference",
      baseline_study = "S1", times = 1
    ),
    "Marginal RMST currently requires piecewise-exponential survival baselines. M-spline fits still support survival, risk, and time-specific marginal hazard measures.",
    fixed = TRUE
  )
})

test_that("M-spline marginal RMST ratios are rejected", {
  fit <- make_mspline_marginal_survival_fit()
  expect_error(
    marginal_effects(
      fit, data.frame(x = 0), measure = "rmst_ratio",
      baseline_study = "S1", times = 1
    ),
    "Marginal RMST currently requires piecewise-exponential survival baselines. M-spline fits still support survival, risk, and time-specific marginal hazard measures.",
    fixed = TRUE
  )
})

test_that("M-spline RMST rejection precedes posterior extraction", {
  fit <- make_mspline_marginal_survival_fit()
  fit$fit$draws <- function(...) stop("posterior draws should not be requested")
  expect_error(
    cpaic:::.cpaic_marginal_survival_predictions(
      fit, X = matrix(0, nrow = 1L), weights = 1,
      treatments = c("Placebo", "A"), baseline_study = "S1",
      times = 1, type = "rmst"
    ),
    "Marginal RMST currently requires piecewise-exponential survival baselines. M-spline fits still support survival, risk, and time-specific marginal hazard measures.",
    fixed = TRUE
  )
})

test_that("M-spline fits retain pointwise marginal survival measures", {
  fit <- make_mspline_marginal_survival_fit()
  survival <- marginal_effects(
    fit, data.frame(x = 0), measure = "survival_difference",
    baseline_study = "S1", times = 1
  )
  risk <- marginal_effects(
    fit, data.frame(x = 0), measure = "risk_difference",
    baseline_study = "S1", times = 1
  )
  hazard <- marginal_effects(
    fit, data.frame(x = 0), measure = "time_specific_hazard_ratio",
    baseline_study = "S1", times = 1
  )
  expect_s3_class(survival, "cpaic_effects")
  expect_s3_class(risk, "cpaic_effects")
  expect_s3_class(hazard, "cpaic_effects")
  expect_true(all(is.finite(c(survival$estimate, risk$estimate, hazard$estimate))))
})
