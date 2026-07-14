test_that("cmlnmr recovers component effects when identified", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("randtoolbox")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "cmdstan not installed")

  set.seed(7)
  beta <- c(A = 0.5, B = 0.4, C = 0.3)
  gamma <- c(A = 0.2, B = 0, C = 0.1)
  Cmat <- build_C_matrix(c("Placebo", "A", "A+B", "C", "A+B+C"),
                         inactive = "Placebo")
  gen <- function(study, trt, n, mux1, mu0) {
    x1 <- rnorm(n, mux1, 1)
    tc <- Cmat[trt, ]
    eta <- mu0 + 0.3 * x1 + sum(tc * beta) + sum(tc * gamma) * x1
    data.frame(.study = study, .trt = trt, .y = rbinom(n, 1, plogis(eta)),
               x1 = x1)
  }
  ipd <- rbind(gen("S1", "Placebo", 400, 0, -0.2), gen("S1", "A", 400, 0, -0.2),
               gen("S2", "A", 400, 0.3, -0.1), gen("S2", "A+B", 400, 0.3, -0.1),
               gen("S5", "Placebo", 400, 0.1, 0), gen("S5", "C", 400, 0.1, 0))
  mk <- function(study, trt, n, mux1, mu0) {
    d <- gen(study, trt, n, mux1, mu0)
    data.frame(.study = study, .trt = trt, r = sum(d$.y), n = nrow(d),
               x1_mean = mean(d$x1), x1_sd = sd(d$x1))
  }
  agd <- rbind(mk("S3", "Placebo", 500, -0.3, 0), mk("S3", "C", 500, -0.3, 0),
               mk("S4", "A+B", 500, 0.5, 0.1), mk("S4", "A+B+C", 500, 0.5, 0.1))

  fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
                chains = 2, iter_warmup = 400, iter_sampling = 400, seed = 1,
                prior_reg_sd = 2.5)
  ce <- component_effects(fit)
  est <- setNames(ce$estimate, ce$component)

  # Each truth inside the 95% credible interval, point estimates in range.
  for (cc in c("A", "B", "C")) {
    row <- ce[ce$component == cc, ]
    expect_true(beta[[cc]] >= row$lower && beta[[cc]] <= row$upper,
                info = paste("component", cc))
  }
  expect_lt(max(fit$fit$summary("beta")$rhat), 1.1)

  # relative_effects() must work on the Bayesian fit (no $bridge slot), and
  # must be told which population to report in.
  expect_error(relative_effects(fit), "Specify `newdata`")
  re <- relative_effects(fit, newdata = data.frame(x1 = 0))
  expect_s3_class(re, "cpaic_effects")
  expect_true(all(c("treatment", "estimate", "lower", "upper") %in% names(re)))
  expect_true(nrow(re) >= 1)
  expect_error(additivity_test(fit), "LOO")

  # The effect is population-specific: gamma must actually be used, so a
  # different target population must give a different answer. Component A has
  # a positive interaction (gamma_A = 0.2), so its effect must increase with
  # x1. Reporting C %*% beta alone would make these identical.
  re0 <- relative_effects(fit, newdata = data.frame(x1 = 0), backtransf = FALSE)
  re1 <- relative_effects(fit, newdata = data.frame(x1 = 1), backtransf = FALSE)
  a0 <- re0$estimate[re0$treatment == "A"]
  a1 <- re1$estimate[re1$treatment == "A"]
  expect_false(isTRUE(all.equal(a0, a1)))
  expect_gt(a1, a0)
  # theta_A(x) = beta_A + gamma_A * x, so the increment recovers gamma_A.
  expect_equal(a1 - a0, 0.2, tolerance = 0.15)

  # Component effects in a target population differ from the main effects.
  ce1 <- component_effects(fit, newdata = data.frame(x1 = 1))
  expect_false(isTRUE(all.equal(ce$estimate, ce1$estimate)))
})

test_that("cmlnmr runs for gaussian, poisson, and survival families", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("randtoolbox")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "cmdstan not installed")

  Cmat <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")
  set.seed(3)

  # gaussian
  genN <- function(study, trt, n, mu0) {
    x1 <- rnorm(n); tc <- Cmat[trt, ]
    y <- mu0 + 0.3 * x1 + sum(tc * c(A = 0.5, B = 0.4)) + rnorm(n)
    data.frame(.study = study, .trt = trt, .y = y, x1 = x1)
  }
  ipdN <- rbind(genN("S1", "Placebo", 200, 0), genN("S1", "A", 200, 0),
                genN("S2", "A", 200, 0), genN("S2", "A+B", 200, 0))
  agdN <- data.frame(.study = "S3", .trt = c("Placebo", "A"),
                     .y = c(0.1, 0.6), se = c(0.1, 0.1),
                     x1_mean = c(0, 0), x1_sd = c(1, 1))
  fN <- cmlnmr(ipdN, agdN, effect_modifiers = "x1", inactive = "Placebo",
               family = "gaussian", chains = 1, iter_warmup = 250,
               iter_sampling = 250, seed = 1)
  expect_true(all(is.finite(component_effects(fN)$estimate)))

  # poisson
  genP <- function(study, trt, n, mu0) {
    x1 <- rnorm(n); tc <- Cmat[trt, ]
    lam <- exp(mu0 - 1 + 0.2 * x1 + sum(tc * c(A = 0.4, B = 0.3)))
    data.frame(.study = study, .trt = trt, .y = rpois(n, lam),
               .exposure = 1, x1 = x1)
  }
  ipdP <- rbind(genP("S1", "Placebo", 300, 0), genP("S1", "A", 300, 0))
  agdP <- data.frame(.study = "S2", .trt = c("A", "A+B"),
                     r = c(120, 150), E = c(300, 300),
                     x1_mean = c(0, 0), x1_sd = c(1, 1))
  fP <- cmlnmr(ipdP, agdP, effect_modifiers = "x1", inactive = "Placebo",
               family = "poisson", chains = 1, iter_warmup = 250,
               iter_sampling = 250, seed = 1)
  expect_true(all(is.finite(component_effects(fP)$estimate)))

  # survival (exponential). The exact likelihood needs individual event and
  # censoring times, so an aggregate survival arm is supplied as RECONSTRUCTED
  # pseudo-individual rows (as multinma::set_agd_surv expects from a digitized
  # Kaplan-Meier curve). Event counts plus person-time cannot identify the
  # individual likelihood, which is precisely what made the old approximation
  # biased upward.
  genS <- function(study, trt, n, mu0) {
    x1 <- rnorm(n); tc <- Cmat[trt, ]
    rate <- exp(mu0 - 1 + 0.2 * x1 + sum(tc * c(A = 0.4, B = 0.3)))
    t <- rexp(n, rate)
    data.frame(.study = study, .trt = trt, .y = 1L, .time = t, x1 = x1)
  }
  ipdS <- rbind(genS("S1", "Placebo", 300, 0), genS("S1", "A", 300, 0))
  recS <- rbind(genS("S2", "A", 250, 0), genS("S2", "A+B", 250, 0))
  agdS <- recS
  agdS$x1_mean <- 0
  agdS$x1_sd <- 1
  fS <- cmlnmr(ipdS, agdS, effect_modifiers = "x1", inactive = "Placebo",
               family = "survival", chains = 1, iter_warmup = 250,
               iter_sampling = 250, seed = 1, n_int = 16)
  expect_true(all(is.finite(component_effects(fS)$estimate)))
})

test_that("cmlnmr fits a piecewise-exponential (flexible) survival baseline", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("randtoolbox")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "cmdstan not installed")

  set.seed(5)
  Cmat <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")
  beta <- c(A = 0.5, B = 0.3)
  gen <- function(study, trt, n, mux1) {
    x1 <- rnorm(n, mux1, 1)
    tc <- Cmat[trt, ]
    lp <- 0.3 * x1 + sum(tc * beta)
    t <- (rexp(n) / (0.03 * exp(lp)))^(1 / 1.3)   # Weibull baseline
    data.frame(.study = study, .trt = trt, .y = as.integer(t <= 24),
               .time = pmin(t, 24), x1 = x1)
  }
  ipd <- rbind(gen("S1", "Placebo", 400, 0), gen("S1", "A", 400, 0),
               gen("S2", "A", 400, 0.3), gen("S2", "A+B", 400, 0.3))
  cuts <- c(6, 12)
  # Reconstructed pseudo-individual survival rows, as multinma::set_agd_surv
  # expects from a digitized Kaplan-Meier curve. Event counts plus person-time
  # cannot identify the exact individual likelihood; that approximation was
  # biased upward by about 36%.
  mk <- function(study, trt, n, mux1) {
    d <- gen(study, trt, n, mux1)
    d$x1_mean <- 0
    d$x1_sd <- 1
    d
  }
  agd <- rbind(mk("S3", "Placebo", 600, -0.2), mk("S3", "A+B", 600, 0.4))

  fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
                family = "survival", cut_points = cuts, chains = 2,
                iter_warmup = 400, iter_sampling = 400, seed = 1, n_int = 16)
  ce <- component_effects(fit)
  expect_true(all(is.finite(ce$estimate)))
  expect_true(all(ce$estimate > 0))                 # both log-HRs positive
  expect_lt(max(fit$fit$summary("beta")$rhat), 1.2)
  # Combined effect A+B in a sensible region of the truth (0.8).
  ab <- sum(ce$estimate)
  expect_gt(ab, 0.4)
  expect_lt(ab, 1.2)
})

test_that("Gaussian copula correlates the integration points", {
  skip_if_not_installed("randtoolbox")
  R <- matrix(c(1, 0.8, 0.8, 1), 2)
  Xc <- .cpaic_integration_points(c(0, 0), c(1, 1), 4096, cor = R)
  Xi <- .cpaic_integration_points(c(0, 0), c(1, 1), 4096, cor = NULL)
  expect_gt(stats::cor(Xc[, 1], Xc[, 2]), 0.6)        # copula induces it
  expect_lt(abs(stats::cor(Xi[, 1], Xi[, 2])), 0.1)   # independent otherwise
  expect_equal(colMeans(Xc), c(0, 0), tolerance = 0.05)  # margins preserved
  expect_equal(apply(Xc, 2, stats::sd), c(1, 1), tolerance = 0.05)
})

test_that("cmlnmr fits an M-spline survival baseline", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("randtoolbox")
  skip_if_not_installed("splines2")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "cmdstan not installed")

  set.seed(5)
  Cmat <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")
  beta <- c(A = 0.5, B = 0.3)
  gen <- function(study, trt, n, mux1) {
    x1 <- rnorm(n, mux1, 1)
    tc <- Cmat[trt, ]
    lp <- 0.3 * x1 + sum(tc * beta)
    t <- (rexp(n) / (0.03 * exp(lp)))^(1 / 1.3)
    data.frame(.study = study, .trt = trt, .y = as.integer(t <= 24),
               .time = pmin(t, 24), x1 = x1)
  }
  ipd <- rbind(gen("S1", "Placebo", 400, 0), gen("S1", "A", 400, 0),
               gen("S2", "A", 400, 0.3), gen("S2", "A+B", 400, 0.3))
  cuts <- c(4, 8, 12, 16)
  # Reconstructed pseudo-individual survival rows, as multinma::set_agd_surv
  # expects from a digitized Kaplan-Meier curve. Event counts plus person-time
  # cannot identify the exact individual likelihood; that approximation was
  # biased upward by about 36%.
  mk <- function(study, trt, n, mux1) {
    d <- gen(study, trt, n, mux1)
    d$x1_mean <- 0
    d$x1_sd <- 1
    d
  }
  agd <- rbind(mk("S3", "Placebo", 600, -0.2), mk("S3", "A+B", 600, 0.4))

  fit <- cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
                family = "survival", cut_points = cuts, baseline = "mspline",
                n_basis = 4, chains = 2, iter_warmup = 400,
                iter_sampling = 400, seed = 1, n_int = 16)
  ce <- component_effects(fit)
  expect_true(all(is.finite(ce$estimate)))
  expect_true(all(ce$estimate > 0))
  expect_lt(max(fit$fit$summary("beta")$rhat), 1.2)
})

test_that("integration margins respect the covariate's support", {
  skip_if_not_installed("randtoolbox")
  # A binary effect modifier with prevalence 0.3. A normal margin would place
  # a third of the integration points outside {0, 1}, integrating the model
  # over a population that cannot exist.
  pts_norm <- cpaic:::.cpaic_integration_points(
    means = 0.3, sds = sqrt(0.21), n_int = 64, margins = "normal")
  expect_true(mean(pts_norm < 0 | pts_norm > 1) > 0.2)

  pts_bern <- cpaic:::.cpaic_integration_points(
    means = 0.3, sds = NA_real_, n_int = 64, margins = "bernoulli")
  expect_true(all(pts_bern %in% c(0, 1)))
  expect_equal(mean(pts_bern), 0.3, tolerance = 0.05)
})

test_that("binary effect modifiers get a Bernoulli margin by default", {
  ipd <- data.frame(.study = "S1", .trt = "A",
                    x1 = rep(c(0, 1), 50), x2 = rnorm(100))
  m <- cpaic:::.cpaic_guess_margins(ipd, c("x1", "x2"))
  expect_equal(unname(m[["x1"]]), "bernoulli")
  expect_equal(unname(m[["x2"]]), "normal")
})

test_that("copula correlation is pooled within studies, not across them", {
  # Two studies, each with ZERO within-study correlation, but shifted means.
  # Pooling all rows manufactures a large spurious correlation.
  set.seed(3)
  n <- 400
  s1 <- data.frame(.study = "S1", x1 = rnorm(n, -2), x2 = rnorm(n, -2))
  s2 <- data.frame(.study = "S2", x1 = rnorm(n, 2), x2 = rnorm(n, 2))
  ipd <- rbind(s1, s2)

  pooled <- stats::cor(ipd$x1, ipd$x2)
  expect_gt(pooled, 0.5)                       # the artefact

  R <- cpaic:::.cpaic_copula_cor(ipd, c("x1", "x2"), ".study")
  expect_lt(abs(R[1, 2]), 0.15)                # the truth
})

test_that("a supplied correlation matrix must really be a correlation matrix", {
  ipd <- data.frame(.study = "S1", x1 = rnorm(50), x2 = rnorm(50))
  expect_error(
    cpaic:::.cpaic_copula_cor(ipd, c("x1", "x2"), ".study",
                              given = diag(2) * 2),
    "unit diagonal")
  bad <- matrix(c(1, 1.5, 1.5, 1), 2)
  expect_error(
    cpaic:::.cpaic_copula_cor(ipd, c("x1", "x2"), ".study", given = bad),
    "positive definite")
})
