#------------------------------------------------------------------------------
## Test Weir
#------------------------------------------------------------------------------


test_that("flow_depth_weir correctly ", {
  # Bollrich S. 421
  expect_equal(flow_depth_weir(B=20, Q=94.9, w = 6.58, mu = 0.797)$h,
               1.58, tolerance = 0.01)
  expect_equal(flow_depth_weir(B=20, Q=94.9, w = 6.58, mu = 0.797)$v,
               0.58, tolerance = 0.01)
}
)


test_that("flow_weir correctly ", {
  # Bollrich S. 421
  expect_equal(flow_weir(B=20, h=1.58, w = 6.58, mu = 0.797)$Q,
               94.9, tolerance = 0.01)
  expect_equal(flow_weir(B=20, h=1.58, w = 6.58, mu = 0.797)$v,
               0.58, tolerance = 0.01)
}
)
