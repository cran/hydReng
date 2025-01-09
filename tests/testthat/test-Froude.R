#------------------------------------------------------------------------------
## Test Froude
#------------------------------------------------------------------------------

#source: https://www.cedengineering.com/userfiles/Hydraulic%20Jumps%20R1.pdf
test_that("Froude in csArbitrary", {
   x <- c(0,1.8,4.8,6.6)/3.281
   z <- c(0.6,0,0,0.6)/3.281
   csA<-CSarbitrary(x=x, z=z, xb_l=1.8/3.281, xb_r = 4.8/3.281,kSt_B = 40)
   v<-5.55/3.281
   h<-0.6/3.281


   froude_number(csA,v=v,h=h)^2

  # Test expected value with a tolerance of 0.01
  expect_equal(froude_number(csA,v=v,h=h)^2, 2.19, tolerance = 0.01)
})


#source Hager S.139 example 6.2
test_that("Froude in cscircle", {

  csC<-CScircle(Di=0.9)
  h<-0.519
  Q<-0.8
  A<-wetted_area(csC,h=h)
  v<-Q/A

  # Test expected value with a tolerance of 0.01
  expect_equal( froude_number(csC,v=v,h=h), 1, tolerance = 0.01)
})
