#------------------------------------------------------------------------------
# Functions for calculating pressflow
#------------------------------------------------------------------------------

# Lambda ks
# --------------------------------------------------------------
lambda_ks <- function(ks, Di, Re) {

  # Function for searching root
  fcn <- function(lam) {
    return(1/sqrt(lam) + (2*log10((2.51/(Re*sqrt(lam))) +
                                    ((ks/Di)/3.71))))
  }

  # Search root
  lam <- uniroot(fcn, interval = c(0.001, 0.1), check.conv = TRUE)$root

  return(lam)
}


# Lambda kst
# --------------------------------------------------------------
lambda_kst <- function(kst, dhy) {
  return(124.58 / (kst^2) * dhy^(1/3))
}


# Reynolds Zahl
# --------------------------------------------------------------
reynolds <- function(v, Di, nu = 1.14e-6) {
  return((v * Di) / nu)
}


# Equivalent Diameter
# --------------------------------------------------------------
#' @title Equivalent Hydraulic Diameter
#' @name d_aequiv
#' @description Calculates the equivalent hydraulic diameter of a rectangular
#' cross-section given its width and height.
#' @usage d_aequiv(b, h)
#' @param b Width of the rectangle [m].
#' @param h Height of the rectangle [m].
#' @return The equivalent hydraulic diameter [m].
#' @examples
#' d_aequiv(b = 2, h = 1)
#' @export

d_aequiv <- function(b, h) {
  # Calculate hydraulic radius
  rhy <- (b * h) / (2 * (b + h))

  # Convert hydraulic radius to hydraulic diameter
  daeq <- 4 * rhy

  return(daeq)
}


# Continuous friction loss
# --------------------------------------------------------------
xi_r <- function(lambda, L, Di) {
  xi <- lambda * L / Di
  return(xi)
}


# Exit loss (Borda-Carnot, in: Hager S. 36)
# --------------------------------------------------------------
xi_a <- function(A1, A2) {
  xi <- (1 - (A1/A2))^2
  return(xi)
}


# Pressflow velocity to air (case 1)
# --------------------------------------------------------------
pressflow_velocity <- function(z0, z1, h0, Di, b = NULL, h = NULL, L, ks, kst, xi_e = 0.5,
                            nu = 1.14e-6, calc_lam = 'kst') {

  if (!is.null(b)) {
    Di <- d_aequiv(b, h)
  }

  fcn <- function(v) {
    Re <- reynolds(v, Di, nu)

    lam <- if (calc_lam == 'kst') {
      lambda_kst(kst, Di)
    } else if (calc_lam == 'ks') {
      lambda_ks(ks, Di, Re)
    }

    xi <- xi_r(lam, L, Di)

    return(v - sqrt(((z0 - z1 + h0) * (2 * 9.81)) / (1 + xi_e + xi)))
  }

  v <- uniroot(fcn, interval = c(0.001, sqrt(2 * 9.81 * (z0 - z1 + h0))),
               check.conv = TRUE)$root

  return(v)
}



# Pressflow to air (case 1)
# --------------------------------------------------------------
#' @title Flow Under Pressure (Bernoulli)
#' @name pressflow
#' @description Calculates the flow in a pipe or a rectangle under pressure
#'   (Bernoulli). The outlet is not submerged, e.g., the exit loss equals 0.
#' @usage pressflow(z0, z1, h0, Di=NULL, h = NULL, b = NULL, L, ks=NULL, kst,
#'   xi_e = 0.5, nu = 1.14e-6, calc_lam = "kst")
#' @param z0 Absolute height of upper gate – upstream of the inlet [m.a.s.l].
#' @param z1 Absolute height of the pipe/rectangle vertical middle axis at
#'   lower gate [m.a.s.l].
#' @param h0 Water depth upstream of the gate – upstream of the inlet [m].
#' @param Di Diameter of pipe [m]. If Di is specified, h and b must be NULL.
#' @param h Height of rectangle [m]. If h is specified, Di must be NULL.
#' @param b Width of rectangle [m]. If b is specified, Di must be NULL.
#' @param L Length of pipe [m].
#' @param ks Equivalent sand roughness [m].
#' @param kst Roughness [m^(1/3)/s].
#' @param calc_lam Defines if lambda should be calculated with ks or kst.
#' @param xi_e Entrance loss [-]. Default = 0.5.
#' @param nu Kinematic viscosity [m2/s]. Default = 1.14e-6.
#'
#' @return Pressflow returns the flow under pressure:
#' \describe{
#'   \item{Q}{Discharge [m^3/s].}
#'   \item{v}{Flow velocity [m/s].}
#' }
#'
#' @examples
#' # Calculate flow in a pipe under pressure with ks value
#' pressflow(z0 = 415, z1 = 413, h0 = 3, L = 20, Di = 1, ks = 0.01,
#'   calc_lam = "ks")
#'
#' # Calculate flow in rectangle under pressure with kst value
#' pressflow(z0 = 415, z1 = 413, h0 = 3, L = 20, b = 2, h = 1, kst = 60,
#'   calc_lam = "kst")
#'
#' @export




pressflow <- function(z0, z1, h0, Di = NULL, h = NULL, b = NULL, L, ks = NULL,
                            kst, xi_e = 0.5, nu = 1.14e-6, calc_lam = 'kst') {

  if (is.null(b) & !is.null(h)) {
    warning("if h is not NULL, b must not be NULL")
  }

  if (is.null(h) & !is.null(b)) {
    warning("if b is not NULL, h must not be NULL")
  }

  if (any(!is.null(h), !is.null(b)) & !is.null(Di)) {
    warning("if Di is not NULL, h and b must be NULL")
  }

  if (any(is.null(h), is.null(b)) & is.null(Di)) {
    warning("if Di is NULL, h and b must not be NULL")
  }

  if (!is.null(b)) {
    Di <- d_aequiv(b, h)
  }

  v <- pressflow_velocity(z0 = z0, z1 = z1, h0 = h0, Di = Di, h = h, b = b, L = L,
                       kst = kst, ks = ks, xi_e = xi_e, nu = nu, calc_lam = calc_lam)

  q <- if (!is.null(b)) {
    v * h * b
  } else {
    v * (Di / 2)^2 * pi
  }

  return(list(Q = q, v = v))
}


# Backwater Height Upstream A Inlet Under Pressure (Case 1)
# --------------------------------------------------------------
#' @title Backwater Height Upstream A Inlet Under Pressure (Bernoulli)
#' @name pressflow_depth
#' @description Calculates the backwater height upstream of an inlet (pipe or
#' rectangle) under pressure (Bernoulli). The outlet is not submerged, e.g.,
#' the exit loss equals 0.
#' @usage pressflow_depth(
#'   z0, z1, Q, Di = NULL, h = NULL, b = NULL, L, ks = NULL, kst,
#'   xi_e = 0.5, nu = 1.14e-6, calc_lam = "kst"
#' )
#' @param z0 Absolute height of upper gate – upstream of the inlet [m.a.s.l].
#' @param z1 Absolute height of the pipe/rectangle vertical middle axis at
#' lower gate [m.a.s.l].
#' @param Q Flow [m^3/s].
#' @param Di Diameter of pipe [m]. If Di is specified, h and b must be NULL.
#' @param h Height of rectangle [m]. If h is specified, Di must be NULL.
#' @param b Width of rectangle [m]. If b is specified, Di must be NULL.
#' @param L Length of pipe [m].
#' @param ks Equivalent sand roughness [m].
#' @param kst Roughness [m^(1/3)/s].
#' @param calc_lam Defines if lambda should be calculated with ks or kst.
#' @param xi_e Entrance loss [-]. Default = 0.5.
#' @param nu Kinematic viscosity [m^2/s]. Default = 1.14e-6.
#'
#' @return Returns the backwater height upstream the inlet:
#' \describe{
#'   \item{h0}{Water depth upstream the inlet [m].}
#'   \item{v}{Flow velocity [m/s].}
#' }
#'
#' @examples
#' # Flow in a pipe under pressure with ks value
#' pressflow_depth(z0 = 415, z1 = 413, Q = 5.18, L = 20, Di = 1,
#'                 ks = 0.01, calc_lam = "ks")
#'
#' # Flow in a rectangle under pressure with kst value
#' pressflow_depth(z0 = 415, z1 = 413, Q = 13.7, L = 20, b = 2, h = 1,
#'                 kst = 60, calc_lam = "kst")
#'
#' @export

pressflow_depth <- function(z0, z1, Q, Di = NULL, h = NULL, b = NULL, L, ks = NULL,
                            kst, xi_e = 0.5, nu = 1.14e-6, calc_lam = 'kst') {

  if (is.null(b) & !is.null(h)) {
    warning("if h is not NULL, b must not be NULL")
  }

  if (is.null(h) & !is.null(b)) {
    warning("if b is not NULL, h must not be NULL")
  }

  if (any(!is.null(h), !is.null(b)) & !is.null(Di)) {
    warning("if Di is not NULL, h and b must be NULL")
  }

  if (any(is.null(h), is.null(b)) & is.null(Di)) {
    warning("if Di is NULL, h and b must not be NULL")
  }

  if (!is.null(b)) {
    Di <- d_aequiv(b, h)
  }

  # Calculate velocity
  v <- if (!is.null(b)) {
    Q / (h * b)
  } else {
    Q / ((Di / 2)^2 * pi)
  }

  lam <- if (calc_lam == 'kst') {
    lambda_kst(kst, Di)
  } else if (calc_lam == 'ks') {
    Re <- reynolds(v, Di, nu)
    lambda_ks(ks, Di, Re)
  }

  xi <- xi_r(lam, L, Di)

  h0 <- (z1 + (v^2 / (2 * 9.81)) * (1 + xi_e + xi)) - z0

  return(list(h0 = h0, v = v))
}



# Backwater Height Upstream A Inlet Under Pressure (Case 2+3)
# --------------------------------------------------------------
#' @title Backwater Height Upstream A Inlet Under Pressure (Bernoulli)
#' @name pressflow_depth_sub
#' @description Calculates the backwater height upstream of an inlet (pipe or
#' rectangle) under pressure (Bernoulli). The outlet is submerged; hence, an
#' exit loss (xi_a) has to be specified.
#' @usage pressflow_depth_sub(
#'   z0, z1, Q, h1, v1, Di = NULL, h = NULL, b = NULL, L, ks = NULL, kst, xi_a,
#'   xi_e = 0.5, nu = 1.14e-6, calc_lam = "kst"
#' )
#' @param z0 Absolute height of upper gate – upstream of the inlet [m.a.s.l].
#' @param z1 Absolute height of the pipe/rectangle vertical middle axis at
#' lower gate [m.a.s.l].
#' @param Q Flow [m^3/s].
#' @param h1 	Water depth downstream the outlet [m].
#' @param v1 Velocity downstream the outlet [m/s].
#' @param Di Diameter of pipe [m]. If Di is specified, h and b must be NULL.
#' @param h Height of rectangle [m]. If h is specified, Di must be NULL.
#' @param b Width of rectangle [m]. If b is specified, Di must be NULL.
#' @param L Length of pipe [m].
#' @param ks Equivalent sand roughness [m].
#' @param kst Roughness [m^(1/3)/s].
#' @param calc_lam Defines if lambda should be calculated with ks or kst.
#' @param xi_a Exit loss, according to Borda-Carnot formula (1 - A1/A2)^2 [-].
#' @param xi_e Entrance loss [-]. Default = 0.5.
#' @param nu Kinematic viscosity [m^2/s]. Default = 1.14e-6.
#'
#' @return Returns the backwater height upstream the inlet:
#' \describe{
#'   \item{h0}{Water depth upstream the inlet [m].}
#'   \item{v}{Flow velocity [m/s].}
#' }
#'
#' @examples
#' # Flow in a pipe under pressure with ks value
#' pressflow_depth_sub(z0=415,z1=413,Q=5.18,h1=2,v1=4,L=20,Di=1,ks=0.01,
#' calc_lam="ks",xi_a=0.5)
#'
#' # Flow in a rectangle under pressure with kst value
#' pressflow_depth_sub(z0=415,z1=413,Q=13.7,h1=2,v1=4,L=20,b=2,h=1,kst=60,
#' calc_lam="kst",xi_a=0.5)
#'
#' @export

pressflow_depth_sub <- function(z0, z1, Q, h1, v1, Di = NULL, h = NULL, b = NULL, L,
                            ks = NULL, kst, xi_a, xi_e = 0.5, nu = 1.14e-6,
                            calc_lam = 'kst') {

  if (is.null(b) & !is.null(h)) {
    warning("if h is not NULL, b must not be NULL")
  }

  if (is.null(h) & !is.null(b)) {
    warning("if b is not NULL, h must not be NULL")
  }

  if (any(!is.null(h), !is.null(b)) & !is.null(Di)) {
    warning("if Di is not NULL, h and b must be NULL")
  }

  if (any(is.null(h), is.null(b)) & is.null(Di)) {
    warning("if Di is NULL, h and b must not be NULL")
  }

  if (!is.null(b)) {
    Di <- d_aequiv(b, h)
  }

  # Calculate velocity
  v <- if (!is.null(b)) {
    Q / (h * b)
  } else {
    Q / ((Di / 2)^2 * pi)
  }

  lam <- if (calc_lam == 'kst') {
    lambda_kst(kst, Di)
  } else if (calc_lam == 'ks') {
    Re <- reynolds(v, Di, nu)
    lambda_ks(ks, Di, Re)
  }

  xi <- xi_r(lam, L, Di)

  h0 <- z1 - z0 + h1 + v1^2 / (2 * 9.81) + v^2 / (2 * 9.81) * (xi_e + xi + xi_a)

  return(list(h0 = h0, v = v))
}
