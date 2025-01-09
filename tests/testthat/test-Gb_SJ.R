#------------------------------------------------------------------------------
## Test Bedload transport according to Smart and Jaeggi
#------------------------------------------------------------------------------

test_that("bedload_SJ calculates correctly", {
  B <- 27
  dm <- 0.11
  d90 <- 0.23
  d30 <- 3200000 * d90 / 4084101  # Calculate d30 based on d90 and constants
  J <- 0.009
  um <- 2.6
  Rs <- 0.8

  # Example from "BAFU(2014):Abschaetzung der mittleren jährlichen
  #Geschiebelieferung in Vorfluter, Praxishilfe", S. 57
  result <- bedload_SJ(d30 = d30, dm = dm, d90 = d90, J = J, Rs = Rs,
                  um = um, B = B, t_crit = 0.023)

  # Test expected value with a tolerance of 0.01
  expect_equal(result, 84.8576, tolerance = 0.01)
})
