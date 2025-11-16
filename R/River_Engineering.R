#------------------------------------------------------------------------------
# Block size according to Stevens and Simons (Bezzola, 2012)
#------------------------------------------------------------------------------

#' Calculate dimensions of rip rap block size
#'
#' Calculates the dimensions and mass of a rip rap block based on slope
#' geometry, water table levels, and material properties.
#'
#' @usage block_size(h, h_z, J, gamma, psi, geo = NULL, S = 1.15, Theta_c = 0.047,
#'   s = 2.65, ret = "all")
#'
#' @param h Numeric. Global maximum water table level above riverbed [m].
#' @param h_z Numeric. Local water table level above the regarded block [m].
#' @param J Numeric. Bottom slope [-].
#' @param gamma Numeric or NULL. Angle of bank slope [degrees]. Use NULL if
#'   specifying \code{geo}.
#' @param geo Numeric vector of length 2 or NULL. Slope geometry as a triangle:
#'   c(vertical length, horizontal length) [-]. If given, \code{gamma} is ignored.
#' @param psi Numeric. Inner friction angle [degrees]. Values between 50 and 55
#'   are recommended (Bezzola 2012).
#' @param S Numeric. Safety factor, default is 1.15 [-].
#' @param Theta_c Numeric. Shear stress parameter, default is 0.047 [-].
#' @param s Numeric. Relative density of blocks, default is 2.65 [-].
#' @param ret Character. Result to return: \code{"all"} (default), \code{"D"}, or \code{"b"}.
#'
#'
#' @return
#' If \code{ret = "all"}, returns a list with:
#' \item{D}{Diameter of block [m]}
#' \item{m}{Mass of block [kg]}
#' \item{a}{a-axis length [m]}
#' \item{b}{b-axis length [m]}
#' \item{c}{c-axis length [m]}
#'
#' Otherwise returns the requested single value:
#' \itemize{
#'  \item \code{"D"} Diameter of block [m]
#'  \item \code{"b"} b-axis length [m]
#' }
#'
#' @references
#' Bezzola (2012). Flussbau, Vorlesungsmanuskript, ETH Zuerich
#'
#' @examples
#' # Calculate block size at bottom of slope with given slope angle
#' block_size(h = 5, h_z = 5, J = 0.0015, gamma = 33.69, psi = 50)
#'
#' # Calculate block size with slope geometries 2:3
#' block_size(h = 5, h_z = 5, J = 0.0015, gamma = NULL, psi = 50, geo = c(2, 3))
#'
#' # Calculate block size at middle of slope with slope geometries 2:3
#' block_size(h = 5, h_z = 2.5, J = 0.0015, gamma = NULL, psi = 50, geo = c(2, 3))
#' @export


block_size <- function(h, h_z, J, gamma, psi,
                       geo = NULL, S = 1.15, Theta_c = 0.047,
                       s = 2.65, ret = "all") {

  if (!is.null(geo)) {
    if (!is.null(gamma)) {
      warning("If 'geo' is specified, set gamma = NULL")
    }

    xrad <- atan(geo[1] / geo[2])
    gamma <- 180 * xrad / pi
  }

  if (h_z > 0.77 * h) {
    H <- 0.77 * h
  } else {
    H <- h_z
  }

  D <- (H * J) /
    (Theta_c * (s - 1) * cos(gamma * pi / 180) *
       (1 / S - S * tan(gamma * pi / 180)^2 / tan(psi * pi / 180)^2))

  m <- round((1 / 6) * D^3 * pi * s * 1e3)

  if (ret == "all") {
    res <- list(
      D = round(D, 2),
      m = m,
      a = round(D / 0.68, 2),
      b = round(D / 0.91, 2),
      c = round(D / 1.30, 2)
    )
    return(res)

  } else if (ret == "D") {
    return(round(D, 2))

  } else if (ret == "b") {
    return(D / 0.91)
  }
}

#------------------------------------------------------------------------------
# Filterlayer
#------------------------------------------------------------------------------
#' Calculate grain size distribution of a filter layer
#'
#' Tool to calculate the range of the grain size distribution of a filter layer.
#'
#' @param d15B Numeric. d15 of block [mm].
#' @param d50B Numeric. d50 of block [mm].
#' @param d15U Numeric. d15 of soil [mm].
#' @param d50U Numeric. d50 of soil [mm].
#' @param d85U Numeric. d85 of soil [mm].
#' @param dmax Numeric. Maximum grain diameter of filter layer [mm].
#' @param plot Logical. If TRUE, the results are plotted (default is TRUE).
#' @param fuller Logical. If TRUE, adds curves of Fuller distributions with
#'   exponents 0.5 < q < 1.5 to the plot. For an ideal grain size distribution,
#'   q is estimated as 0.5 (default is FALSE).
#'
#'
#' @return A list with the following components:
#'   \item{d15min}{Minimum d15 of filter layer [mm].}
#'   \item{d15max}{Maximum d15 of filter layer [mm].}
#'   \item{d50min}{Minimum d50 of filter layer [mm].}
#'   \item{d50max}{Maximum d50 of filter layer [mm].}
#'   \item{d85min}{Minimum d85 of filter layer [mm].}
#'
#' @importFrom graphics abline text
#' @importFrom grDevices rainbow terrain.colors
#'
#' @examples
#' # Calculate range of the grain size distribution
#' filterlayer(1000, 1500, 5, 10, 20, 400)
#'
#' # Calculate range of the grain size distribution and add Fuller curves
#' filterlayer(1000, 1500, 5, 10, 20, 400, fuller = TRUE)
#'
#' @export
#'

filterlayer <- function(d15B, d50B, d15U, d50U, d85U, dmax = 400, plot = TRUE,
                        fuller = FALSE) {
  d15F_min1 <- d15B / 20
  d15F_max1 <- d15B / 5
  d85F_min1 <- d15B / 5
  d50F_min1 <- d50B / 25

  d15F_max2_1 <- d15U * 20
  d15F_min2 <- d15U * 5
  d15F_max2_2 <- d85U * 5
  d50F_max2 <- d50U * 25

  d15F_max2 <- min(d15F_max2_1, d15F_max2_2)

  if (plot) {
    plot(
      c(d15F_min1, d15F_max1, d85F_min1, d50F_min1),
      c(0.15, 0.15, 0.85, 0.5),
      xlim = c(0.01, 1500),
      ylim = c(0.001, 1),
      log = "x",
      col = "blue",
      xlab = "grain size [mm]",
      ylab = "%",
      main = "Specification of filterlayer",
      xaxt = "n"
    )
    axis(1, at = 10 ^ (-4:3), labels = 10 ^ (-4:3))
    points(c(d15F_min2, d15F_max2, d50F_max2), c(0.15, 0.15, 0.5),
           col = "blue", pch = 16)
    points(c(d15B, d50B), c(0.15, 0.5), pch = 15)
    points(c(d15U, d50U, d85U), c(0.15, 0.5, 0.85), pch = 17)
    abline(v = c(1:10 %o% 10 ^ (-3:3)), col = "gray")
    abline(h = c(0.2, 0.4, 0.6, 0.8), col = "gray")

    if (fuller) {
      vfull <- seq(0.5, 1.5, 0.2)
      for (i in seq_along(vfull)) {
        vfuller <- numeric(dmax)
        for (j in seq_len(dmax)) {
          vfuller[j] <- (j / dmax) ^ vfull[i]
        }
        lines(1:dmax, vfuller, col = rainbow(length(vfull))[i])
      }
    }

    polygon(
      c(max(d15F_min1, d15F_min2), min(d15F_max2, d15F_max1),
        min(d15F_max2, d15F_max1), max(d15F_min1, d15F_min2)),
      c(0.16, 0.16, 0.14, 0.14),
      border = "red"
    )
    polygon(c(d50F_min1, d50F_max2, d50F_max2, d50F_min1),
            c(0.51, 0.51, 0.49, 0.49), border = "red")
    polygon(c(d85F_min1, 1500, 1500, d85F_min1),
            c(0.86, 0.86, 0.84, 0.84), border = "red")

    polygon(c(0.002, 0.002, 0.06, 0.06), c(0, 1, 1, 0),
            col = terrain.colors(4, alpha = 0.2)[1])
    polygon(c(0.06, 0.06, 2, 2), c(0, 1, 1, 0),
            col = terrain.colors(4, alpha = 0.2)[2])
    polygon(c(2, 2, 60, 60), c(0, 1, 1, 0),
            col = terrain.colors(4, alpha = 0.2)[3])
    polygon(c(60, 60, 10000, 10000), c(0, 1, 1, 0),
            col = terrain.colors(4, alpha = 0.2)[4])

    text(0.001, 0.95, "Clay")
    text(0.01, 0.95, "Silt")
    text(0.4, 0.95, "Sand")
    text(8, 0.95, "Gravel")
    text(400, 0.95, "Cobble")

    if (!fuller) {
      legend(
        "bottomleft",
        legend = c("Soil", "Block", "Criterion 1", "Criterion 2",
                   "Distribution Range"),
        pch = c(17, 15, 1, 16, 0),
        col = c("black", "black", "blue", "blue", "red"),
        cex = 0.8
      )
    } else {
      legend(
        "bottomleft",
        legend = c(
          "Soil", "Block", "Criterion 1", "Criterion 2", "Distribution Range",
          "Fuller Exponent", "0.5", "0.7", "0.9", "1.1", "1.3", "1.5"
        ),
        pch = c(17, 15, 1, 16, 0, rep(NA, 7)),
        lty = c(rep(NA, 6), rep(1, 6)),
        col = c("black", "black", "blue", "blue", "red", NA, rainbow(6)),
        cex = 0.8
      )
    }
  }

  list(
    d15min = d15F_min1,
    d15max = d15F_max2,
    d50min = d50F_min1,
    d50max = d50F_max2,
    d85min = d85F_min1
  )
}



#------------------------------------------------------------------------------
# Superelevation of water table in curve (Bezzola 2012, Kap. 11.4)
#------------------------------------------------------------------------------
#' Superelevation of water table in curve
#'
#' Calculates the superelevation of the water table in a river curve.
#'
#' @param w Numeric. Horizontal sole width [m].
#' @param rm Numeric. Curve radius from center to the middle of the river [m].
#' @param v Numeric. Flow velocity [m/s].
#' @param S Numeric. Safety factor, default is 1.5.
#'
#' @return Numeric. The difference between mean water level and superelevation [m].
#'
#' @references
#' Bezzola (2012). Flussbau, Vorlesungsmanuskript, ETH Zuerich
#'
#' @examples
#' # Calculate superelevation
#' wt_sup(w = 30, rm = 200, v = 5)
#' @export

wt_sup <- function(w, rm, v, S = 1.5) {
  dh <- S * (w / rm) * (v^2 / (2 * 9.81))
  return(dh)
}

#------------------------------------------------------------------------------
# Shear stress calculation
#------------------------------------------------------------------------------

#' Shear stress, shear velocity, and dimensionless shear stress
#'
#' Calculates shear stress, shear velocity, and dimensionless shear stress
#' based on water depth, slope, and grain size.
#'
#' @param h0 Numeric. Total water depth [m].
#' @param J Numeric. Bottom slope [-].
#' @param dm Numeric or NULL. Median grain size (\code{d50}) of sediment [mm].
#' @param h Numeric or NULL. Local water depth at the point of interest [m].
#'   If \code{NULL}, considered equal to \code{h0}.
#' @param rho Numeric. Density of water [kg/m3], default is 1000.
#'
#' @return A named list with components:
#' \item{tau}{Shear stress [N/m2].}
#' \item{U}{Shear velocity [m/s].}
#' \item{tau_st}{Dimensionless shear stress [-], if \code{dm} is provided,
#' otherwise \code{NA}.}
#'
#' @references
#' Bezzola (2012). Flussbau, Vorlesungsmanuskript, ETH Zuerich
#'
#' @examples
#' # Calculate shear stress at bank bottom
#' shear_str(h0 = 3.31, J = 0.0022)$tau
#'
#' # Calculate shear stress at bank middle
#' shear_str(h0 = 3.31, J = 0.0022, h = 1.6)$tau
#'
#' # Calculate dimensionless shear stress
#' shear_str(h0 = 3.31, J = 0.0022, dm = 100)$tau_st
#' @export


shear_str <- function(h0, J, dm = NULL, h = NULL, rho = 1000) {

  if (is.numeric(h)) {
    tau <- rho * 9.81 * h0 * J * (1 - ((h0 - h) / h0))
  } else {
    tau <- rho * 9.81 * h0 * J
  }

  # Shear velocity
  U <- sqrt(tau / rho)

  # Dimensionless shear stress
  if (is.numeric(dm)) {
    tau_st <- tau / (1650 * 9.81 * (dm / 1000))
    return(list(tau = tau, U = U, tau_st = tau_st))
  } else {
    return(list(tau = tau, U = U))
  }
}
