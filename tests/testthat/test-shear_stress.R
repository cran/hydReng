#------------------------------------------------------------------------------
## Test shear_stress
#------------------------------------------------------------------------------

test_that("shear_stress calculates correctly", {
  result <- shear_str(h0 = 3.31, J = 0.0022)
  expect_equal(result$tau, 71.43, tolerance = 0.01)
  expect_equal(result$U, 0.267, tolerance = 0.01)

})

