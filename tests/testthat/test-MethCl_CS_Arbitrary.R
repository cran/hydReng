#------------------------------------------------------------------------------
## Test CSarbitrary Methods
#------------------------------------------------------------------------------

# Setup function
setUp <- function() {
  # Create profiles
  x <- c(0, 4, 9, 13)
  z <- c(2, 0, 0, 2)
  csAr_Bollrich3_6_1_6 <<- CSarbitrary(x = x, z = z, xb_l = 4, xb_r = 9,
                                       kSt_B = 35, kSt_l = 45, kSt_r = 45)

  x <- seq(-sqrt(20), sqrt(20), length.out = 41)
  z <- 0.1 * x^2 - 2
  csAr_Bollrich3_6_1_7 <<- CSarbitrary(x = x, z = z, kSt_B = 40)

  x <- c(-0.85, 3, 15, 18.85)
  z <- c(3.85, 0, 0, 3.85)

  # Initialize the CSarbitrary object
  csAr_KOHS_Schaechen <<- CSarbitrary(
    x = x, z = z, xb_l = 3, xb_r = 15,
    kSt_B = 45
  )
}


# Call setup function
setUp()

#------------------------------------------------------------------------------
## Test access and geometry
#------------------------------------------------------------------------------

## Test access
test_that("Test Arbitrary Access", {
  expect_equal(kSt_B(csAr_Bollrich3_6_1_7), 40)
})

## Test replacements that are not allowed
test_that("Test Arbitrary Replacement", {
  expect_error(x(csAr_Bollrich3_6_1_7) <- 1:3)
})

## Test unique x values are not allowed
test_that("Test unique x values", {
  x <- c(0, 0, 2, 4)
  z <- c(5, 0, 0, 5)
  expect_error(CSarbitrary(x = x, z = z, kSt_B = 40))
})

#------------------------------------------------------------------------------
## Hydraulics
#------------------------------------------------------------------------------
test_that("wetted_area for Bollrich3_6_1_7 calculates correctly", {
  expect_equal(wetted_area(csAr_Bollrich3_6_1_7, h = 2), 11.92, tolerance = 0.001)
})

test_that("wetted_perimeter for Bollrich3_6_1_7 calculates correctly", {
  expect_equal(wetted_perimeter(csAr_Bollrich3_6_1_7, h = 2), 10.02,
               tolerance = 0.001)
})


test_that("wetted_area for Bollrich3_6_1_6 calculates correctly", {
  expect_equal(wetted_area(csAr_Bollrich3_6_1_6, h = 2), 18, tolerance = 0.001)
})


test_that("flow_velocity for Bollrich3_6_1_6 calculates correctly", {
  expect_equal(flow_velocity(csAr_Bollrich3_6_1_6, h = 2, J = 0.0001,
                             method = "Einstein"), 0.482, tolerance = 0.001)
})

test_that("flow_depth for Bollrich3_6_1_6 calculates correctly", {
  expect_equal(flow_depth(csAr_Bollrich3_6_1_6, Q = 8.677, J = 0.0001,
                                  method = "Einstein", ret = "h"),
               2, tolerance = 0.001)
  expect_equal(flow_depth(csAr_Bollrich3_6_1_6, Q = 8.677, J = 0.0001,
                                  method = "Einstein")$v, 0.482, tolerance = 0.001)
  expect_equal(flow_depth(csAr_Bollrich3_6_1_6, Q = 8.677, J = 0.0001,
                                  method = "Einstein")$A, 18, tolerance = 0.001)
  expect_equal(flow_depth(csAr_Bollrich3_6_1_6, Q = 8.677, J = 0.0001,
                                  method = "Einstein")$P, 13.94, tolerance = 0.001)
  expect_equal(flow_depth(csAr_Bollrich3_6_1_6, Q = 8.677, J = 0.0001,
                                  method = "Einstein")$kSt_m, 40.66306, tolerance = 0.001)
})

test_that("mean_roughness for Bollrich3_6_1_6 calculates correctly", {
  expect_equal(mean_roughness(csAr_Bollrich3_6_1_6, h=2),
               40.66306, tolerance = 0.001)
})

test_that("flow for Bollrich3_6_1_6 calculates correctly", {
  expect_equal(flow(csAr_Bollrich3_6_1_6, h=2, J=0.0001, method="Einstein",
                    ret="Q"),
               8.677, tolerance = 0.001)
  expect_equal(flow(csAr_Bollrich3_6_1_6, h=2, J=0.0001,
                                      method="Einstein")$v,0.482, tolerance = 0.001)
  expect_equal(flow(csAr_Bollrich3_6_1_6, h=2, J=0.0001,
                                      method="Einstein")$A,18, tolerance = 0.001)
  expect_equal(flow(csAr_Bollrich3_6_1_6, h=2, J=0.0001,
                                      method="Einstein")$kSt_m,40.66306, tolerance = 0.001)
})


test_that("flow_max for Bollrich3_6_1_6 calculates correctly", {
  expect_equal(flow_max(csAr_Bollrich3_6_1_6, J=0.0001, method="Einstein",
                                      ret="Qmax"),
               8.677, tolerance = 0.001)
  expect_equal(flow_max(csAr_Bollrich3_6_1_6, J=0.0001, method="Einstein",
                                 ret="v"),
               0.482, tolerance = 0.001)
  expect_equal(flow_max(csAr_Bollrich3_6_1_6, J=0.0001, method="Einstein",
                                 ret="hmax"),
               2, tolerance = 0.001)
  expect_equal(flow_max(csAr_Bollrich3_6_1_6, J=0.0001,
                                 method="Einstein")$A,18, tolerance = 0.001)

})



test_that("flow_max_freeboard calculations are accurate", {
  result_flow_max_freeboard <- flow_max_freeboard(
    csAr_KOHS_Schaechen,
    J = 2.2e-2,
    type = "KOHS",
    sigma_wz = 0,
    fv = TRUE,
    ft = 0.5,
    fe = NULL,
    method = "Strickler"
  )
  expect_equal(result_flow_max_freeboard$fe, 2.57, tolerance = 0.01,
               label = "flow_max_freeboard 'fe' value is within tolerance")
  expect_equal(result_flow_max_freeboard$v, 7.1, tolerance = 0.01,
               label = "flow_max_freeboard 'v' value is within tolerance")
  expect_equal(result_flow_max_freeboard$Qmax, 120, tolerance = 0.01,
               label = "flow_max_freeboard 'Qmax' value is within tolerance")
  expect_equal(result_flow_max_freeboard$hmax, 1.27, tolerance = 0.01,
               label = "flow_max_freeboard 'hmax' value is within tolerance")


})
