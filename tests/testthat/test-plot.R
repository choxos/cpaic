# Plots. We do not snapshot images (no vdiffr dependency); we assert that each
# function returns a ggplot object and that ggplot2 can actually build it, which
# catches missing aesthetics, unbound variables, and bad scales.

skip_if_no_ggplot <- function() skip_if_not_installed("ggplot2")

expect_ggplot <- function(p) {
  expect_s3_class(p, "gg")
  expect_no_error(ggplot2::ggplot_build(p))
  invisible(p)
}

bin_net <- function(ipd = TRUE) {
  if (ipd) {
    cpaic_network(cpaic_bin_agd, ipd = cpaic_bin_ipd, sm = "OR",
                  family = "binomial", ipd_covariates = "x1",
                  inactive = "Placebo")
  } else {
    cpaic_network(cpaic_bin_agd, sm = "OR", inactive = "Placebo")
  }
}

# Network plot -----------------------------------------------------------------

test_that("plot.cpaic_network returns a buildable ggplot", {
  skip_if_no_ggplot()
  expect_ggplot(plot(bin_net()))
  expect_ggplot(plot(bin_net(ipd = FALSE)))
  expect_ggplot(plot(bin_net(), weight_edges = FALSE, show_bridges = FALSE))
  expect_ggplot(plot(bin_net(), nudge = 0))
})

test_that("the network plot shows the disconnection and names the bridges", {
  skip_if_no_ggplot()
  net <- bin_net()
  conn <- cpaic_connectivity(net)
  expect_false(conn$connected)
  p <- plot(net)
  # Nodes are colored by sub-network, so every sub-network must appear.
  node_layer <- ggplot2::ggplot_build(p)$plot$layers[[2]]$data
  expect_setequal(unique(node_layer$subnetwork),
                  paste("Sub-network", seq_len(conn$n_subnetworks)))
  # The bridging components are named for the reader.
  expect_true(grepl("bridging components", p$labels$subtitle, fixed = TRUE))
  for (cc in conn$bridging_components) {
    expect_true(grepl(cc, p$labels$subtitle, fixed = TRUE))
  }
})

test_that("plot.cpaic_network validates its arguments", {
  skip_if_no_ggplot()
  net <- bin_net()
  expect_error(plot(net, weight_edges = "yes"), "TRUE or FALSE")
  expect_error(plot(net, show_bridges = NA), "TRUE or FALSE")
  expect_error(plot(net, nudge = c(1, 2)), "single finite number")
})

# Forest plots ------------------------------------------------------------------

test_that("forest() builds for bridges, effect tables, and components", {
  skip_if_no_ggplot()
  br <- cnma_bridge(bin_net(ipd = FALSE))
  expect_ggplot(forest(br))
  expect_ggplot(forest(br, what = "component"))
  expect_ggplot(forest(relative_effects(br)))
  expect_ggplot(forest(component_effects(br)))
  expect_ggplot(forest(br, order = "alphabetical"))
  expect_ggplot(forest(br, order = "none"))
  expect_ggplot(forest(br, ref_line = NA))
  expect_ggplot(forest(br, backtransf = FALSE))
  # And through the plot() generic.
  expect_ggplot(plot(br))
  expect_ggplot(plot(relative_effects(br)))
})

test_that("forest() SHOWS non-estimable contrasts rather than dropping them", {
  skip_if_no_ggplot()
  # A + B and A + B + C versus Placebo need the B component, which no
  # within-study contrast identifies here.
  net <- cpaic_network(cpaic_bin_agd[c(1, 3), ], sm = "OR",
                       inactive = "Placebo")
  br <- suppressWarnings(cnma_bridge(net))
  re <- suppressWarnings(relative_effects(br))
  expect_true(anyNA(re$estimate))

  p <- suppressWarnings(forest(br))
  expect_ggplot(p)
  built <- ggplot2::ggplot_build(p)
  n_na <- sum(is.na(re$estimate))

  # Every contrast keeps its row, estimable or not: none is silently dropped.
  expect_setequal(levels(built$plot$data$.label),
                  paste0(re$treatment, " vs ", re$comparator))
  # A "not estimable" text mark is drawn for each of them ...
  txt <- built$data[[length(built$data)]]
  expect_equal(nrow(txt), n_na)
  # ... and they are named in the caption.
  expect_true(grepl("Not estimable", p$labels$caption))
  for (t in re$treatment[is.na(re$estimate)]) {
    expect_true(grepl(t, p$labels$caption, fixed = TRUE))
  }

  # show_na = FALSE drops them, and then the caption is gone.
  p2 <- suppressWarnings(forest(br, show_na = FALSE))
  expect_ggplot(p2)
  expect_null(p2$labels$caption)
})

test_that("forest() facets an all-contrasts table by comparator", {
  skip_if_no_ggplot()
  br <- cnma_bridge(bin_net(ipd = FALSE))
  re <- relative_effects(br, all_contrasts = TRUE)
  expect_gt(length(unique(re$comparator)), 1L)
  p <- forest(re)
  expect_ggplot(p)
  expect_s3_class(p$facet, "FacetWrap")
  # Faceted rows are labelled by treatment; the comparator is the panel.
  expect_setequal(levels(p$data$.label), unique(re$treatment))
})

test_that("forest() rejects a bad ref_line", {
  skip_if_no_ggplot()
  br <- cnma_bridge(bin_net(ipd = FALSE))
  expect_error(forest(br, ref_line = c(0, 1)), "single number")
  expect_error(forest(br, ref_line = "zero"), "single number")
  # 0 is not the null value of an odds ratio, and a log axis cannot show it.
  expect_error(forest(br, ref_line = 0), "must be\\s+positive")
  expect_ggplot(forest(br, ref_line = 0, backtransf = FALSE))
})

test_that("forest() draws ratio measures on a log axis", {
  skip_if_no_ggplot()
  br <- cnma_bridge(bin_net(ipd = FALSE))
  # OR is back-transformed, so the x scale is log10 ...
  expect_true(any(vapply(forest(br)$scales$scales,
                         function(s) identical(s$trans$name, "log-10") ||
                           identical(s$transform$name, "log-10"),
                         logical(1))))
  # ... but the link scale is linear.
  expect_false(any(vapply(forest(br, backtransf = FALSE)$scales$scales,
                          function(s) identical(s$trans$name, "log-10") ||
                            identical(s$transform$name, "log-10"),
                          logical(1))))
})

# Edge influence -----------------------------------------------------------------

test_that("plot_edge_influence() builds and flags zero-influence IPD edges", {
  skip_if_no_ggplot()
  br <- cnma_bridge(bin_net())
  p <- plot_edge_influence(br, treatment = "A+B+C")
  expect_ggplot(p)
  expect_true(all(c("AgD", "IPD") %in%
                    levels(ggplot2::ggplot_build(p)$plot$data$edge_type)))

  # S4 (A+B+D vs A+B) carries IPD but cannot inform the D-free contrast
  # A vs Placebo, so it must be flagged rather than passed over.
  p2 <- suppressWarnings(plot_edge_influence(br, treatment = "A"))
  expect_ggplot(p2)
  dat <- ggplot2::ggplot_build(p2)$plot$data
  dead <- dat$studlab[dat$zero_influence]
  expect_true(length(dead) > 0)
  expect_true(grepl("Zero influence", p2$labels$caption))
})

# Bayesian plots -----------------------------------------------------------------

skip_if_no_stan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("randtoolbox")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "cmdstan not installed")
}

# One small fit, reused by every Bayesian plot test in this file.
mlnmr_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(11)
    Cmat <- build_C_matrix(c("Placebo", "A", "A+B", "A+B+C"),
                           inactive = "Placebo")
    beta <- c(A = 0.5, B = 0.4, C = 0.3)
    gam <- c(A = 0.3, B = 0, C = -0.5)
    gen <- function(study, trt, n, mux1, mu0) {
      x1 <- stats::rnorm(n, mux1, 1)
      tc <- Cmat[trt, ]
      eta <- mu0 + 0.3 * x1 + sum(tc * beta) + sum(tc * gam) * x1
      data.frame(.study = study, .trt = trt,
                 .y = stats::rbinom(n, 1, stats::plogis(eta)), x1 = x1)
    }
    mk <- function(study, trt, n, mux1, mu0) {
      d <- gen(study, trt, n, mux1, mu0)
      data.frame(.study = study, .trt = trt, r = sum(d$.y), n = nrow(d),
                 x1_mean = mean(d$x1), x1_sd = stats::sd(d$x1))
    }
    ipd <- rbind(gen("S1", "Placebo", 200, 0, -0.2), gen("S1", "A", 200, 0, -0.2),
                 gen("S2", "A", 200, 0.4, 0), gen("S2", "A+B", 200, 0.4, 0),
                 gen("S5", "A+B", 200, -0.2, 0.1),
                 gen("S5", "A+B+C", 200, -0.2, 0.1))
    agd <- rbind(mk("S3", "Placebo", 300, -0.3, 0), mk("S3", "A", 300, -0.3, 0),
                 mk("S4", "A+B", 300, 0.5, 0.1),
                 mk("S4", "A+B+C", 300, 0.5, 0.1))
    cache <<- suppressWarnings(
      cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
             chains = 2, iter_warmup = 200, iter_sampling = 200, seed = 3))
    cache
  }
})

test_that("rank_probs() returns a proper rank distribution", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  rp <- rank_probs(fit, newdata = data.frame(x1 = 0.3), what = "component")
  expect_s3_class(rp, "cpaic_rank_probs")
  # Every element's rank probabilities sum to one.
  tot <- tapply(rp$probability, rp$element, sum)
  expect_equal(unname(as.numeric(tot)), rep(1, length(tot)), tolerance = 1e-8)

  cum <- rank_probs(fit, newdata = data.frame(x1 = 0.3), what = "component",
                    cumulative = TRUE)
  expect_true(all(cum$probability >= -1e-8 & cum$probability <= 1 + 1e-8))
  top <- cum$probability[cum$rank_position == max(cum$rank_position)]
  expect_equal(top, rep(1, length(top)), tolerance = 1e-8)
  expect_error(rank_probs(fit, newdata = data.frame(x1 = 0), cumulative = NA),
               "TRUE or FALSE")
})

test_that("rankogram and cumulative rank plots build", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  nd <- data.frame(x1 = 0.3)
  expect_ggplot(plot(rank_probs(fit, newdata = nd, what = "component")))
  expect_ggplot(plot(rank_probs(fit, newdata = nd, what = "component",
                                cumulative = TRUE)))
  expect_ggplot(plot(rank_probs(fit, newdata = nd)))
  expect_ggplot(plot(cpaic_ranks(fit, newdata = nd, what = "component")))
  expect_ggplot(plot(cpaic_ranks(fit, newdata = nd), metric = "mean_rank"))
})

test_that("plot_rank_curve() traces the hierarchy across target populations", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  vals <- seq(-1, 1, by = 0.5)
  p <- plot_rank_curve(fit, em = "x1", values = vals, what = "component")
  expect_ggplot(p)
  # One point per (component, target value): the hierarchy really is a curve.
  expect_setequal(unique(ggplot2::ggplot_build(p)$plot$data$target), vals)

  # It also accepts the data frame that rank_curve() returns.
  rc <- rank_curve(fit, em = "x1", values = vals, what = "component")
  expect_ggplot(plot_rank_curve(rc))
  expect_ggplot(plot_rank_curve(fit, em = "x1", values = vals,
                                metric = "p_best"))
  expect_error(plot_rank_curve(fit), "`em` and `values` are required")
  expect_error(plot_rank_curve(data.frame(a = 1)), "cmlnmr\\(\\) fit")
})

test_that("plot_estimability() maps estimability over target populations", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  p <- plot_estimability(fit, em = "x1", values = seq(-1, 1, by = 0.5))
  expect_ggplot(p)
  lv <- levels(ggplot2::ggplot_build(p)$plot$data$identified_by)
  expect_setequal(lv, c("IPD (within-study)", "Aggregate (ecological)",
                        "Not estimable"))
  expect_error(plot_estimability(fit, em = "nope", values = 0),
               "effect modifiers")
  expect_error(plot_estimability(fit, em = "x1", values = numeric(0)),
               "finite target values")
})

test_that("DIC, leverage, prior-posterior, and integration plots build", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  expect_ggplot(plot(dic(fit)))
  expect_ggplot(plot(dic(fit), dic(fit), labels = c("A", "B")))
  expect_ggplot(plot_leverage(fit))
  expect_ggplot(plot_leverage(fit, dic_contours = c(1, 3)))
  expect_ggplot(plot_prior_posterior(fit))
  expect_ggplot(plot_prior_posterior(fit, prior = "gamma"))
  expect_ggplot(plot_integration_error(fit))
  expect_ggplot(plot_integration_error(fit, int_thin = 16))

  expect_error(plot(dic(fit), "not a dic"), "cpaic_dic")
  expect_error(plot_prior_posterior(fit, prior = "banana"),
               "character vector with elements")
  # tau exists only under random effects.
  expect_error(plot_prior_posterior(fit, prior = "tau"),
               "character vector with elements")
})

test_that("plot.cpaic_mlnmr gives MCMC diagnostics via bayesplot", {
  skip_if_no_stan()
  skip_if_not_installed("bayesplot")
  fit <- mlnmr_fixture()
  expect_ggplot(plot(fit, type = "trace"))
  expect_ggplot(plot(fit, type = "density"))
  expect_ggplot(plot(fit, type = "hist"))
  expect_ggplot(plot(fit, type = "rhat"))
  expect_ggplot(plot(fit, type = "neff"))
  # mcmc_pairs() returns a bayesplot grid, not a ggplot.
  expect_s3_class(plot(fit, type = "pairs", pars = c("beta[1]", "beta[2]")),
                  "bayesplot_grid")
  expect_error(plot(fit, type = "trace", pars = character(0)),
               "non-empty character vector")
})

test_that("forest() works on a cmlnmr fit in a named target population", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  expect_ggplot(forest(fit, newdata = data.frame(x1 = 0.3)))
  expect_ggplot(forest(fit, what = "component",
                       newdata = data.frame(x1 = 0.3)))
  # Relative effects under population adjustment have no population-free value.
  expect_error(forest(fit), "newdata")
})

test_that("survival curves build, with and without the KM overlay", {
  skip_if_no_stan()
  set.seed(4)
  Cmat <- build_C_matrix(c("Placebo", "A", "A+B"), inactive = "Placebo")
  sg <- function(study, trt, n, mux1) {
    x1 <- stats::rnorm(n, mux1, 1)
    tc <- Cmat[trt, ]
    eta <- 0.2 * x1 + sum(tc * c(A = 0.4, B = 0.3))
    tt <- stats::rexp(n, rate = 0.15 * exp(eta))
    cens <- stats::runif(n, 1, 12)
    data.frame(.study = study, .trt = trt, .time = pmin(tt, cens),
               .y = as.integer(tt <= cens), x1 = x1)
  }
  ipd <- rbind(sg("S1", "Placebo", 120, 0), sg("S1", "A", 120, 0))
  agd <- rbind(sg("S2", "A", 100, 0.5), sg("S2", "A+B", 100, 0.5))
  agd$x1_mean <- stats::ave(agd$x1, agd$.trt, FUN = mean)
  agd$x1_sd <- stats::ave(agd$x1, agd$.trt, FUN = stats::sd)
  agd$x1 <- NULL

  fit <- suppressWarnings(
    cmlnmr(ipd, agd, effect_modifiers = "x1", inactive = "Placebo",
           family = "survival", cut_points = c(3),
           chains = 2, iter_warmup = 150, iter_sampling = 150, seed = 9,
           n_int = 16))
  expect_ggplot(plot_survival(fit, ndraws = 25))
  expect_ggplot(plot_survival(fit, ndraws = 25) + geom_km(fit))
  expect_error(plot_survival(fit, level = 2), "in \\(0, 1\\)")
  expect_error(plot_survival(fit, times = c(-1, 1)), "positive and finite")

  # A leverage plot needs a saturated model, which censored data do not have.
  expect_error(plot_leverage(fit), "saturated model")
  # And the aggregate survival contribution is not an integrated mean outcome.
  expect_error(plot_integration_error(fit), "not supported for survival")
})

test_that("survival helpers reject non-survival fits", {
  skip_if_no_stan()
  fit <- mlnmr_fixture()
  expect_error(plot_survival(fit), "family = \"survival\"")
  expect_error(geom_km(fit), "family = \"survival\"")
})
