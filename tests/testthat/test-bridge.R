test_that("build_C_matrix codes components and inactive correctly", {
  C <- build_C_matrix(c("A", "B", "A + B", "placebo"), inactive = "placebo")
  expect_equal(colnames(C), c("A", "B"))
  expect_equal(unname(C["A + B", ]), c(1L, 1L))
  expect_equal(unname(C["placebo", ]), c(0L, 0L))
})

test_that("connectivity flags a connected, identifiable network", {
  skip_if_not_installed("netmeta")
  data("Linde2016", package = "netmeta")
  net <- cpaic_network(Linde2016, treat1 = "treat1", treat2 = "treat2",
                       TE = "lnOR", seTE = "selnOR", studlab = "id",
                       sm = "OR", inactive = "Placebo")
  conn <- cpaic_connectivity(net)
  expect_true(conn$connected)
  expect_equal(conn$n_subnetworks, 1L)
  expect_true(conn$identifiable)
})

test_that("a disconnected network is bridged through shared components", {
  skip_if_not_installed("netmeta")
  # subnet1 {Placebo, A, B}; subnet2 {A+B, A+B+C}; no shared treatment.
  agd <- data.frame(
    studlab = c("S1", "S2", "S3"),
    treat1  = c("A", "B", "A+B+C"),
    treat2  = c("Placebo", "Placebo", "A+B"),
    TE      = c(0.5, 0.3, 0.2),
    seTE    = c(0.1, 0.1, 0.1)
  )
  net <- cpaic_network(agd, sm = "OR", inactive = "Placebo")
  conn <- cpaic_connectivity(net)

  expect_false(conn$connected)
  expect_equal(conn$n_subnetworks, 2L)
  expect_setequal(conn$bridging_components, c("A", "B"))
  expect_true(conn$identifiable)

  br <- cnma_bridge(net, common = TRUE, random = FALSE)
  ce <- component_effects(br)
  est <- setNames(ce$estimate, ce$component)
  expect_equal(unname(est[c("A", "B", "C")]), c(0.5, 0.3, 0.2),
               tolerance = 1e-6)

  # Cross-subnetwork contrast recovered additively.
  te <- br$fit$TE.common
  expect_equal(unname(te["A+B+C", "Placebo"]), 1.0, tolerance = 1e-6)
  expect_equal(unname(te["A+B+C", "A"]), 0.5, tolerance = 1e-6)
})

test_that("non-bridgeable disconnected network is reported as unidentifiable", {
  skip_if_not_installed("netmeta")
  # Two isolated placebo-anchored stars sharing NO components.
  agd <- data.frame(
    studlab = c("S1", "S2"),
    treat1  = c("A", "C"),
    treat2  = c("PlaceboX", "PlaceboY"),
    TE      = c(0.5, 0.4),
    seTE    = c(0.1, 0.1)
  )
  net <- cpaic_network(agd, sm = "OR")
  conn <- cpaic_connectivity(net)
  expect_false(conn$connected)
  expect_false(conn$identifiable)
})
