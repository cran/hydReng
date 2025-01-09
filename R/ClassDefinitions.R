# Define Class Unions
setClassUnion("numeric_or_logical", c("numeric", "logical"))
setClassUnion("numeric_or_null", c("numeric", "NULL"))

# constructor CSarbitrary
#------------------------------------------------------------------------------
#' @title CSarbitrary Class
#' @aliases CSarbitrary-class
#' @description Defines a cross-section class with arbitrary geometry for
#'   hydraulic calculations. For single open channels only, avoid geometries
#'   with multiple channels.
#' @slot x A numeric vector of x-coordinates [m].
#' @slot z A numeric vector of z-coordinates [m].
#' @slot xb_l x-coordinate of the left bank bottom [m].
#' @slot xb_r x-coordinate of the right bank bottom [m].
#' @slot kSt_B Roughness of the channel bed [m^(1/3)/s].
#' @slot kSt_l Roughness of the left bank [m^(1/3)/s].
#' @slot kSt_r Roughness of the right bank [m^(1/3)/s].
#' @importFrom methods setClass new validObject
#' @examples
#' # Define sample cross-section data
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- new("CSarbitrary", x = x, z = z, xb_l = 4, xb_r = 9,
#'           kSt_B = 35, kSt_l = 45, kSt_r = 45)
#' @export

CSarbitrary <- setClass(
  "CSarbitrary",
  slots = c(
    x = "numeric",
    z = "numeric",
    xb_l = "numeric_or_null",
    xb_r = "numeric_or_null",
    kSt_B = "numeric_or_null",
    kSt_l = "numeric_or_null",
    kSt_r = "numeric_or_null"
  ),
  prototype = list(
    kSt_B = NULL,
    xb_l = NULL,
    xb_r = NULL,
    kSt_l = NULL,
    kSt_r = NULL
  )
)

# constructor CScircle
#------------------------------------------------------------------------------

#' @title CScircle Class
#' @aliases CScircle-class
#' @description Defines a cross-section class with circular geometry for
#'   hydraulic calculations.
#' @slot Di Diameter of the pipe [m].
#' @slot kSt Roughness of the pipe according to Strickler [m^(1/3)/s].
#' @slot ks Roughness of the pipe according to Prandtl-Coolebrook-White [mm]
#'  (SIA 190)
#' @examples
#' csC <- CScircle(Di = 1, kSt = 75)
#' csC <- CScircle(Di = 1, ks = 1.5)
#' @export

CScircle <- setClass(
  "CScircle",
  slots = c(
    Di = "numeric",
    kSt = "numeric_or_null",
    ks = "numeric_or_null"
  ),
  prototype = list(
    kSt = NULL,
    ks = NULL
  )
)







# Set Generics
#------------------------------------------------------------------------------
setGeneric("return_valid_object", function(object, value) {
  standardGeneric("return_valid_object")
})

setGeneric("x", function(object, ...) standardGeneric("x"))
setGeneric("x<-", function(object, value) standardGeneric("x<-"))
setGeneric("z", function(object, ...) standardGeneric("z"))
setGeneric("z<-", function(object, value) standardGeneric("z<-"))
setGeneric("xb_l", function(object, ...) standardGeneric("xb_l"))
setGeneric("xb_l<-", function(object, value) standardGeneric("xb_l<-"))
setGeneric("xb_r", function(object, ...) standardGeneric("xb_r"))
setGeneric("xb_r<-", function(object, value) standardGeneric("xb_r<-"))
setGeneric("kSt_B", function(object, ...) standardGeneric("kSt_B"))
setGeneric("kSt_B<-", function(object, value) standardGeneric("kSt_B<-"))
setGeneric("kSt_l", function(object, ...) standardGeneric("kSt_l"))
setGeneric("kSt_l<-", function(object, value) standardGeneric("kSt_l<-"))
setGeneric("kSt_r", function(object, ...) standardGeneric("kSt_r"))
setGeneric("kSt_r<-", function(object, value) standardGeneric("kSt_r<-"))


setGeneric("Di", function(object, ...) standardGeneric("Di"))
setGeneric("Di<-", function(object, value) standardGeneric("Di<-"))
setGeneric("kSt", function(object, ...) standardGeneric("kSt"))
setGeneric("kSt<-", function(object, value) standardGeneric("kSt<-"))
setGeneric("ks", function(object, ...) standardGeneric("ks"))
setGeneric("ks<-", function(object, value) standardGeneric("ks<-"))



setGeneric(
  "wetted_area",
  function(object, h, ret = "A") standardGeneric("wetted_area")
)

setGeneric(
  "wetted_perimeter",
  function(object, h, ret= "P") standardGeneric("wetted_perimeter")
)
setGeneric("mean_roughness", function(object, h)
  standardGeneric("mean_roughness"))

setGeneric("flow_velocity", function(object,h, J, method="Strickler", nu = 1.14e-6, ...)
  standardGeneric("flow_velocity"))

setGeneric("flow_depth", function(object,  Q, J, method="Strickler",
                                          ret="all", plot=FALSE)
  standardGeneric("flow_depth"))

setGeneric("WL_coords", function(object, ...) standardGeneric("WL_coords"))

setGeneric("EL_coords", function(object, ...) standardGeneric("EL_coords"))

setGeneric("froude_number", function(object, v, h)
  standardGeneric("froude_number"))

setGeneric("flow", function(object, h, J, method = "Strickler", ret = "all",
                                              plot = FALSE)
             standardGeneric("flow"))

setGeneric("flow_max", function(object, J, method= "Strickler",
                                         ret = "all", plot = FALSE)
  standardGeneric("flow_max"))

setGeneric("flow_max_freeboard", function(object, J, type = "KOHS", sigma_wz = 0, fw = TRUE, fv = FALSE, ft = 0,
                                                   fe = NULL, fe_min = 0, fe_max = Inf, method = "Strickler",
                                                   ret = "all", plot = FALSE)
  standardGeneric("flow_max_freeboard"))

setGeneric("par_fill", function(object,J, method="Strickler") standardGeneric("par_fill"))

