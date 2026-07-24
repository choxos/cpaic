test_that("package loads and engine dependencies are available", {
  # The hard dependencies: the component-NMA engine, the MAIC weights, and the
  # Sobol' points every non-Gaussian cmlnmr() integration needs. multinma is a
  # Suggests, used only by the correlation parity test, so it is deliberately
  # absent from this list.
  expect_true(requireNamespace("netmeta", quietly = TRUE))
  expect_true(requireNamespace("maicplus", quietly = TRUE))
  expect_true(requireNamespace("randtoolbox", quietly = TRUE))
  expect_true(requireNamespace("igraph", quietly = TRUE))
  expect_true(requireNamespace("loo", quietly = TRUE))
})
