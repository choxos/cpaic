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
