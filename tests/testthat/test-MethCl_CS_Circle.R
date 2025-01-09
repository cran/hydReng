#------------------------------------------------------------------------------
## Test CScircleMethods
#------------------------------------------------------------------------------


# Setup: create test objects
setUp <- function() {
  # Create profiles
  csC_Skript_ZHW_2006 <<- CScircle(Di = 0.6, ks = 1.5, kSt = 85)
  csC_Hager <<- CScircle(Di = 0.7, kSt = 85)
  csC_Bollrich<<-CScircle(Di=1.2,ks = 1.5,kSt=75)
  csC_Bollrich_2<<-CScircle(Di=3,ks = 1.5)

}

# Setup the test objects
setUp()

# --------------------------------------------------------------------------------------
# Test access and geometry
# --------------------------------------------------------------------------------------

test_that("Circle: access properties correctly", {
  csC_Skript <- csC_Skript_ZHW_2006
  expect_equal(Di(csC_Skript), 0.6)
  expect_equal(kSt(csC_Skript), 85)
  expect_equal(ks(csC_Skript), 1.5)
})

test_that("Circle: replacing properties raises exceptions", {
  csC_Skript <- csC_Skript_ZHW_2006
  expect_error(Di(csC_Skript) <- -10)
})

# --------------------------------------------------------------------------------------
# Test hydraulics
# --------------------------------------------------------------------------------------

test_that("Circle: hydraulics using Prandtl-Colebrook-White method", {
  csC_Skript <- csC_Skript_ZHW_2006
  expect_equal(
    flow(csC_Skript, h = 0.356, J = 15e-3, method = "Prandtl-Coolebrook-White", ret = "Q"),
    0.5,
    tolerance = 0.1
  )
  expect_equal(
    flow_depth(csC_Skript, Q = 0.5, J = 15e-3, method = "Prandtl-Coolebrook-White", ret = "h"),
    0.356,
    tolerance = 0.1
  )
})

test_that("Circle: hydraulics using Strickler method", {
  expect_equal(flow_depth(csC_Hager, Q = 0.46, J = 0.004)$h, 0.43, tolerance = 0.02)
  expect_equal(flow(csC_Hager, h = 0.4359357, J = 0.004)$Q, 0.46, tolerance = 0.02)
})


test_that("Circle: hydraulics using Prandtl-Colebrook-White (Bollrich S.260)", {
  expect_equal(wetted_area(csC_Bollrich,h=1.2),1.13,tolerance=0.1)
  expect_equal(wetted_perimeter(csC_Bollrich,h=1.2),3.76,tolerance=0.1)

  expect_equal(flow(csC_Bollrich,J=0.001,h=0.84,method = "Prandtl-Coolebrook-White")$Q,0.94,tolerance=0.1)
  expect_equal(flow_velocity(csC_Bollrich, J=0.001, h=0.84, method = "Prandtl-Coolebrook-White"),1.15,tolerance=0.1)
  expect_equal(flow_depth(csC_Bollrich, J=0.001, Q=0.9958771, method = "Prandtl-Coolebrook-White")$h,0.84,tolerance=0.1)
  expect_equal(flow_max(csC_Bollrich, J=0.001, method = "Prandtl-Coolebrook-White")$Qmax,1.28,tolerance=0.1)

}
)

