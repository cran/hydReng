#------------------------------------------------------------------------------
## Test Gate
#------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------
# Test flow_gate based on Bollrich examples
#--------------------------------------------------------------------------------------

test_that("flow_gate works correctly based on Bollrich examples", {
  # Bollrich S. 385, free flow
  expect_equal(flow_gate(a = 0.4, h0 = 5.6, B = 5.9, alpha = 90, h2 = NULL, ret = "Q"),
               14.84, tolerance = 0.01)
  expect_equal(flow_gate(a = 0.4, h0 = 5.6, B = 5.9, alpha = 90, h2 = 2, ret = "Q"),
               14.84, tolerance = 0.01)

  # Bollrich S. 385, submerged flow
  expect_equal(flow_gate(a = 0.4, h0 = 5.6, B = 5.9, alpha = 90, h2 = 2.8, ret = "Q"),
               11.85, tolerance = 0.01)
})

#--------------------------------------------------------------------------------------
# Test flow_depth_gate based on Bollrich examples
#--------------------------------------------------------------------------------------

test_that("flow_depth_gate works correctly based on Bollrich examples", {
  # Bollrich S. 385, free flow
  expect_equal(flow_depth_gate(a = 0.4, Q = 14.84, B = 5.9, alpha = 90, h2 = NULL, ret = "h0"),
               5.6, tolerance = 0.01)
  expect_equal(flow_depth_gate(a = 0.4, Q = 14.84, B = 5.9, alpha = 90, h2 = 2, ret = "h0"),
               5.6, tolerance = 0.01)

  # Bollrich S. 385, submerged flow
  expect_equal(flow_depth_gate(a = 0.4, Q = 11.85, B = 5.9, alpha = 90, h2 = 2.8, ret = "h0"),
               5.6, tolerance = 0.01)
})

#--------------------------------------------------------------------------------------
# Test mu_gate based on a specific example
#--------------------------------------------------------------------------------------

test_that("mu_gate works correctly for a specific example", {
  a <- 1
  h0 <- a / 0.725
  alpha <- 30

  psi <- psi_gate(a, h0, alpha)
  mu <- mu_gate(psi, a, h0)

  expect_equal(mu, 0.65, tolerance = 0.01)
})
