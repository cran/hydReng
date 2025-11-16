#------------------------------------------------------------------------------
## Test scour
#------------------------------------------------------------------------------

test_that("scour_vert calculates correctly", {
  result <- scour_vert(Q = 4, B = 1, h = 5, h_u = 1.76, d90 = 150, d95 = 200,
                       ful_ov = TRUE)
  expect_equal(result$T0, 7.722174, tolerance = 0.01)
})

test_that("scour_horz calculates correctly", {
  result <- scour_horz(Q = 4, B = 1, h = 5, h_u = 1.76, d90 = 150, a = 1,
                       mu = 0.6)
  expect_equal(result$T0, 10.66, tolerance = 0.01)
})

test_that("scour_curve calculates correctly", {
  A<-135.5
  Fr<-0.52
  h<-3.31
  J<-0.0022
  r<-500
  rm<-530
  d16<-50
  d84<-200
  result<-scour_curve(A=A,Fr=Fr,h=h,J=J,rm=rm,r=r,d16=d16,d84=d84)
  expect_equal(result$T0, 0.3366091, tolerance = 0.01)
})

test_that("scour_groyne calculates correctly", {
  v<-2.7
  Fr<-0.52
  h<-3.31
  J<-0.0022
  L<-5
  d16<-50
  dm<-80
  d84<-200
  Ks<-0.82
  delta<-60

  result<-scour_groyne(v=v,Fr=Fr, h=h,J=J,L=L,d16=d16,dm=dm,d84=d84,Ks=Ks,
                       delta=delta)


  expect_equal(result$T0, 5.168923, tolerance = 0.01)
})

