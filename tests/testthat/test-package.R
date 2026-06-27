test_that("package loads and engine dependencies are available", {
  expect_true(requireNamespace("netmeta", quietly = TRUE))
  expect_true(requireNamespace("maicplus", quietly = TRUE))
  expect_true(requireNamespace("multinma", quietly = TRUE))
})
