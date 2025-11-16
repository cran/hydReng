#------------------------------------------------------------------------------
## Test filterlayer
#------------------------------------------------------------------------------

test_that("filterlayer calculates correctly", {
  result <- filterlayer(1000, 1500, 5, 10, 20, 400)
  expect_equal(result$d15min, 50, tolerance = 0.01)
  expect_equal(result$d85min, 200, tolerance = 0.01)

})

