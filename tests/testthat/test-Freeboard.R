#------------------------------------------------------------------------------
## Test freeboard
#------------------------------------------------------------------------------


# Test for freeboard with ft=0
test_that("freeboard (ft=0) calculations are accurate", {
  expect_equal(
    freeboard(v = 7.1, h = 1.3, sigma_wz = 0, fv = TRUE, ft = 0),
    2.57, tolerance = 0.005,
    label = "freeboard (v=7.1, h=1.3, ft=0) is accurate"
  )
})

# Test for freeboard with ft=0.5
test_that("freeboard (ft=0.5) calculations are accurate", {
  expect_equal(
    freeboard(v = 7.1, h = 1.3, sigma_wz = 0, fv = TRUE, ft = 0.5),
    2.62, tolerance = 0.001,
    label = "freeboard (v=7.1, h=1.3, ft=0.5) is accurate"
  )
})

# Test for freeboard with sigma_wz=0.3
test_that("freeboard with sigma_wz=0.3 calculations are accurate", {
  expect_equal(
    freeboard(v = 4.56, h = 1.36, sigma_wz = 0.3, fv = TRUE, ft = 0),
    1.11, tolerance = 0.001,
    label = "freeboard (v=4.56, h=1.36, sigma_wz=0.3, ft=0) is accurate"
  )
})
