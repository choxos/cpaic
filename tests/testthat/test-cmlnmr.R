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

  # survival (exponential): event indicator + time
  genS <- function(study, trt, n, mu0) {
    x1 <- rnorm(n); tc <- Cmat[trt, ]
    rate <- exp(mu0 - 1 + 0.2 * x1 + sum(tc * c(A = 0.4, B = 0.3)))
    t <- rexp(n, rate)
    data.frame(.study = study, .trt = trt, .y = 1L, .time = t, x1 = x1)
  }
  ipdS <- rbind(genS("S1", "Placebo", 300, 0), genS("S1", "A", 300, 0))
  agdS <- data.frame(.study = "S2", .trt = c("A", "A+B"),
                     r = c(280, 290), E = c(900, 850),
                     x1_mean = c(0, 0), x1_sd = c(1, 1))
  fS <- cmlnmr(ipdS, agdS, effect_modifiers = "x1", inactive = "Placebo",
               family = "survival", chains = 1, iter_warmup = 250,
               iter_sampling = 250, seed = 1)
  expect_true(all(is.finite(component_effects(fS)$estimate)))
})
