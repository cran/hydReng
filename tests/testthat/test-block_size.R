#------------------------------------------------------------------------------
## Test block_size
#------------------------------------------------------------------------------

test_that("block_size calculates correctly", {
  result <- block_size(h = 5, h_z = 5, J = 0.0015, gamma = 33.69, psi = 50)
  expect_equal(result$D, 0.18, tolerance = 0.01)
})

test_that("block_size calculates correctly", {
  result <- block_size(h = 5, h_z = 5, J = 0.0015, gamma = NULL, psi = 50,
                       geo = c(2, 3))
  expect_equal(result$D, 0.18, tolerance = 0.01)
})
