two_stage_gate_fixture <- function() {
  cols <- list(
    studlab = "studlab", treat1 = "treat1", treat2 = "treat2",
    TE = "TE", seTE = "seTE"
  )
  agd <- data.frame(
    studlab = c("S1", "S2"),
    treat1 = c("A", "B"), treat2 = c("C", "C"),
    TE = c(NA_real_, 0.4), seTE = c(NA_real_, 0.15),
    stringsAsFactors = FALSE
  )
  adjusted <- data.frame(
    studlab = "S1", treat1 = "A", treat2 = "C",
    stringsAsFactors = FALSE
  )
  list(cols = cols, agd = agd, adjusted = adjusted)
}

test_that("two-stage bridges stop by default when aggregate edges remain", {
  d <- two_stage_gate_fixture()

  expect_error(
    cpaic:::.cpaic_two_stage_bridge_gate(
      d$agd, d$adjusted, d$cols, "cstc()", "binomial", FALSE
    ),
    "cannot form a decision-grade component bridge.*retained aggregate-only",
    fixed = FALSE
  )
  expect_error(
    cpaic:::.cpaic_two_stage_bridge_gate(
      d$agd, d$adjusted, d$cols, "cmaic()", "binomial", FALSE
    ),
    "cannot form a decision-grade component bridge.*retained aggregate-only",
    fixed = FALSE
  )
})

test_that("experimental opt-in records every gate failure", {
  d <- two_stage_gate_fixture()

  expect_warning(
    validity <- cpaic:::.cpaic_two_stage_bridge_gate(
      d$agd, d$adjusted, d$cols, "cmaic()", "binomial", TRUE
    ),
    "explicitly exploratory"
  )

  expect_false(validity$decision_grade)
  expect_true(validity$experimental_override)
  expect_length(validity$reasons, 2L)
  expect_match(validity$reasons[[1]], "S2: B vs C", fixed = TRUE)
  expect_match(validity$reasons[[2]], "marginal binomial contrasts", fixed = TRUE)
  expect_equal(validity$retained_aggregate_edges$studlab, "S2")
})

test_that("nonlinear cMAIC is gated even when every edge is adjusted", {
  d <- two_stage_gate_fixture()
  agd <- d$agd[1, , drop = FALSE]

  expect_error(
    cpaic:::.cpaic_two_stage_bridge_gate(
      agd, d$adjusted, d$cols, "cmaic()", "binomial", FALSE
    ),
    "marginal binomial contrasts.*not generally additive"
  )
  expect_warning(
    validity <- cpaic:::.cpaic_two_stage_bridge_gate(
      agd, d$adjusted, d$cols, "cmaic()", "binomial", TRUE
    ),
    "marginal binomial contrasts"
  )
  expect_false(validity$decision_grade)
  expect_true(validity$experimental_override)
  expect_equal(nrow(validity$retained_aggregate_edges), 0L)

  expect_silent(
    gaussian <- cpaic:::.cpaic_two_stage_bridge_gate(
      agd, d$adjusted, d$cols, "cmaic()", "gaussian", FALSE
    )
  )
  expect_true(gaussian$decision_grade)
  expect_false(gaussian$experimental_override)
})

test_that("two-stage edge validation rejects partial multi-arm replacement", {
  d <- two_stage_gate_fixture()
  complete <- d$agd[1, , drop = FALSE]

  expect_silent(
    cpaic:::.cpaic_validate_two_stage_edge(
      complete, d$cols, c("A", "C"), "S1", "cstc()"
    )
  )
  expect_error(
    cpaic:::.cpaic_validate_two_stage_edge(
      complete, d$cols, c("A", "B"), "S1", "cstc()"
    ),
    "adjusted contrast cannot replace a different treatment pair"
  )

  multi_arm <- rbind(
    transform(complete, treat1 = "A", treat2 = "C"),
    transform(complete, treat1 = "B", treat2 = "C")
  )
  expect_error(
    cpaic:::.cpaic_validate_two_stage_edge(
      multi_arm, d$cols, c("A", "C"), "S1", "cmaic()"
    ),
    "partial IPD representation of a multi-arm aggregate study"
  )
})

test_that("IPD-only study edges require explicit, recorded permission", {
  agd <- data.frame(
    studlab = "S2", treat1 = "A", treat2 = "B", TE = 0.2, seTE = 0.1,
    stringsAsFactors = FALSE
  )
  ipd <- data.frame(
    .study = rep("S1", 8), .trt = rep(c("A", "B"), each = 4),
    .y = seq_len(8), x = seq(-1, 1, length.out = 8),
    stringsAsFactors = FALSE
  )
  net <- cpaic_network(
    agd, ipd, sm = "MD", family = "gaussian", ipd_covariates = "x"
  )

  expect_error(
    cpaic:::.cpaic_two_stage_plan(net, NULL, "cstc()"),
    "allow_ipd_only_studies = TRUE", fixed = TRUE
  )
  expect_silent(
    plan <- cpaic:::.cpaic_two_stage_plan(
      net, NULL, "cstc()", allow_ipd_only_studies = TRUE
    )
  )
  expect_identical(plan$ipd_only_studies, "S1")
  expect_true(plan$studies[[1]]$ipd_only)
})

test_that("print methods expose gate status, reasons, and retained edges", {
  d <- two_stage_gate_fixture()
  validity <- suppressWarnings(cpaic:::.cpaic_two_stage_bridge_gate(
    d$agd, d$adjusted, d$cols, "cstc()", "binomial", TRUE
  ))
  common <- list(
    bridge_validity = validity,
    network = list(cols = d$cols),
    effect_modifiers = "x", bridge = "bridge result"
  )
  stc_fit <- structure(
    c(common, list(prognostics = "x")), class = "cpaic_stc"
  )
  maic_fit <- structure(
    c(common, list(ess = c(S1 = 10))), class = "cpaic_maic"
  )

  expect_output(print(stc_fit), "EXPERIMENTAL OVERRIDE ACTIVE", fixed = TRUE)
  expect_output(print(stc_fit), "not decision-grade", fixed = TRUE)
  expect_output(print(stc_fit), "S2: B vs C", fixed = TRUE)
  expect_output(print(stc_fit), "Retained aggregate-only edges: S2: B vs C",
                fixed = TRUE)

  expect_output(print(maic_fit), "EXPERIMENTAL OVERRIDE ACTIVE", fixed = TRUE)
  expect_output(print(maic_fit), "Reasons:", fixed = TRUE)
  expect_output(print(maic_fit), "Retained aggregate-only edges: S2: B vs C",
                fixed = TRUE)

  passed <- common
  passed$bridge_validity <- list(
    decision_grade = TRUE, experimental_override = FALSE,
    reasons = character(), retained_aggregate_edges = d$agd[0, ]
  )
  passed_stc <- structure(
    c(passed, list(prognostics = "x")), class = "cpaic_stc"
  )
  expect_output(print(passed_stc), "Two-stage bridge gate: PASSED", fixed = TRUE)
})
