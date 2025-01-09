#------------------------------------------------------------------------------
## Test Pressflow Air
#------------------------------------------------------------------------------

setup <- function() {
  # Create profiles
  z0 <- 0
  z1 <- 3
  h0 <- 20
  b <- 2
  h <- 3
  Di <- d_aequiv(b, h)
  L <- 13
  ks <- 0.0006
  kst <- 100
  xi_e <- 0.5
  A <- b * h

  # Return the values as a list for easy testing
  list(z0 = z0, z1 = z1, h0 = h0, b = b, h = h, Di = Di, L = L, ks = ks, kst = kst, xi_e = xi_e, A = A)
}

# Test case for pressflow_air
test_that("pressflow_air calculations match expected values", {

  # Initialize setup
  values <- setup()

  # Test lambda_ks function
  expect_equal(lambda_ks(values$ks, values$Di, 14.539 * 2.4 / (1.31e-6)), 0.0144, tolerance = 0.001)

  # Test pressflow_velocity function
  expect_equal(pressflow_velocity(z0 = values$z0, z1 = values$z1, h0 = values$h0,
                               b = values$b, h = values$h, L = values$L, ks = values$ks,
                               xi_e = values$xi_e, calc_lam = 'ks'), 14.539, tolerance = 0.001)

  # Test pressflow function
  expect_equal(pressflow(z0 = values$z0, z1 = values$z1, h0 = values$h0,
                               b = values$b, h = values$h, L = values$L, ks = values$ks,
                               xi_e = values$xi_e, calc_lam = 'ks')$Q, 87.23, tolerance = 0.001)

  # Test pressflow_depth function
  expect_equal(pressflow_depth(z0 = values$z0, z1 = values$z1, Q = 87.23,
                               b = values$b, h = values$h, L = values$L, ks = values$ks,
                               xi_e = values$xi_e, calc_lam = 'ks')$h0, 20, tolerance = 0.001)
})

#------------------------------------------------------------------------------
## Test Pressflow Sub
#------------------------------------------------------------------------------



setup <- function() {
  # Create profiles
  z0 <- 800
  z1 <- 800
  b <- 2
  h <- 3
  L <- 13
  ks <- 0.0006
  xi_e <- 0.3 + 2 * 0.1
  xi_a <- 0
  Q <- 87.23
  h1 <- 3
  v1 <- 14.539
  calc_lam <- 'ks'
  nu <- 1.31 * 10^-6

  # Return the values as a list for easy testing
  list(z0 = z0, z1 = z1, b = b, h = h, L = L, ks = ks, xi_e = xi_e, xi_a = xi_a,
       Q = Q, h1 = h1, v1 = v1, calc_lam = calc_lam, nu = nu)
}

# Test case for pressflow_depth_sub function
test_that("pressflow_sub calculations match expected values", {
  # Initialize setup
  values <- setup()

  # Test pressflow_depth_sub function
  expect_equal(pressflow_depth_sub(z0 = values$z0, z1 = values$z1, Q = values$Q,
                               h1 = values$h1, v1 = values$v1, b = values$b,
                               h = values$h, L = values$L, ks = values$ks,
                               xi_e = values$xi_e, xi_a = values$xi_a,
                               calc_lam = values$calc_lam, nu = values$nu)$h0, 20, tolerance = 0.001)

})
