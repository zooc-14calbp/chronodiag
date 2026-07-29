test_that("EPS reference value is reproduced", {
  expect_equal(eps_value(0.23943661971830986, 18), 0.85, tolerance = 1e-12)
})

test_that("SSS is defined at zero correlation", {
  expect_equal(sss_value(0, 5, 20), 0.25, tolerance = 1e-15)
})

test_that("full subset SSS is one", {
  expect_equal(sss_value(0.3, 20, 20), 1, tolerance = 1e-15)
})

test_that("inadmissible equicorrelation is rejected", {
  expect_error(eps_value(-0.2, 10), "admissible")
})

test_that("identity correlation has full effective dimension", {
  expect_equal(effective_dimension(diag(5)), 5, tolerance = 1e-12)
})
