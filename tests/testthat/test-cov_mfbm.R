test_that("variance mfbm works", {
  expect_equal(var_fbm(1, 0.3, sigma = 1), 1)
})

test_that("cross-covariance mfbm works 1", {
  expect_equal(var_fbm(4, 0.3, sigma = 1), cov_ijts_mfbm(4, 4, 0.3, 0.3, 1, 0))
})

test_that("cross-covariance mfbm works 1", {
  expect_equal(0.5, cov_ijts_mfbm(1, 1, 0.3, 0.5, 0.5, 0.02))
})
