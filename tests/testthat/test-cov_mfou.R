test_that("cov_ijk works 1", {
  expect_equal(round(cov_ijk(sin(0 : 314), cos(0 : 314), 0), 2), 0)
})
test_that("cov_ijk works 2", {
  expect_equal(round(cov_ijk(sin(0 : 314), cos(0 : 314), 1), 2), 0.42)
})
test_that("cov_ijk works 3", {
  expect_equal(round(cov_ijk(sin(0 : 314), cos(0 : 314), - 1), 2), - 0.42)
})
test_that("cov_ijk_mfou works", {
  expect_equal(round(cov_ijk_mfou(1, 2, 1, 2, 0.1, 0.2, 0.5, 0.01, 10), 3),
               - 0.001)
})
test_that("cov_ijk_mfou and cov_jik_mfou work", {
  expect_equal(cov_ijk_mfou(1, 2, 1, 2, 0.1, 0.2, 0.5, 0.01, 10),
               cov_jik_mfou(2, 1, 2, 1, 0.2, 0.1, 0.5, - 0.01, 10))
})
test_that("cov_ijk_mfou and cov_ji0_mfou work", {
  expect_equal(cov_ijk_mfou(1, 2, 1, 2, 0.1, 0.2, 0.5, 0.01, 0),
               cov_ij0_mfou(1, 2, 1, 2, 0.1, 0.2, 0.5, 0.01))
})
