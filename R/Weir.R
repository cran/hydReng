#------------------------------------------------------------------------------
# Functions for calculating flow over a Weir
#------------------------------------------------------------------------------

#' @title Flow Depth At Weir Crest
#' @description Calculates the height difference between the upstream water
#'   level and the weir crest.
#' @usage flow_depth_weir(B, Q, w = Inf, mu = 0.73)
#' @param B Width of the weir [m].
#' @param Q Flow rate [m3/s].
#' @param w Height of the weir crest (upstream) [m]. If w = Inf,
#'   the upstream velocity is considered 0.
#' @param mu Discharge coefficient [-]. Default is 0.73.
#' @return A list with the following components:
#' \describe{
#'   \item{h}{Flow depth over the weir [m].}
#'   \item{v}{Flow velocity [m/s].}
#' }
#' @examples
#' flow_depth_weir(B = 3, Q = 5)
#' flow_depth_weir(B = 3, Q = 5, w = 1)
#' @export
flow_depth_weir <- function(B, Q, w = Inf, mu = 0.73) {
  # Function for finding the root
  fcn <- function(h) {
    v <- Q / ((h + w) * B)
    term1 <- (h + v^2 / (2 * 9.81))^(3/2)
    term2 <- (v^2 / (2 * 9.81))^(3/2)
    Q - (2 / 3 * mu * B * sqrt(2 * 9.81) * (term1 - term2))
  }

  # Search root
  h_max <- 0.5 * 3^(2/3) * (Q / (B * sqrt(9.81) * mu))^(2/3) +
    0.5 * 3^(2/3) * (Q / (B * sqrt(9.81) * mu))^(2/3)
  h <- uniroot(fcn, interval = c(1e-5, h_max), check.conv = TRUE)$root
  v <- Q / ((h + w) * B)

  return(list(h = h, v = v))
}


#' @title Flow Over Weir Crest
#' @description Calculates the flow over a weir crest based on upstream
#'   water level.
#' @usage flow_weir(B, h, w = Inf, mu = 0.73)
#' @param B Width of the weir [m].
#' @param h Height difference between the upstream water level and the
#'   weir crest [m].
#' @param w Height of the weir crest (upstream) [m]. If w = Inf, the
#'   upstream velocity is considered 0.
#' @param mu Discharge coefficient [-]. Default is 0.73.
#' @return A list with the following components:
#' \describe{
#'   \item{Q}{Flow over the weir [m3/s].}
#'   \item{v}{Flow velocity [m/s].}
#' }
#' @examples
#' flow_weir(B = 3, h = 1.2)
#' flow_weir(B = 3, h = 1.2, w = 1)
#' @export
#'
flow_weir <- function(B, h, w = Inf, mu = 0.73) {
  # Function for finding the root
  fcn <- function(Q) {
    v <- Q / ((h + w) * B)
    term1 <- (h + v^2 / (2 * 9.81))^(3/2)
    term2 <- (v^2 / (2 * 9.81))^(3/2)
    Q - (2 / 3 * mu * B * sqrt(2 * 9.81) * (term1 - term2))
  }

  # Search root
  Q_max <- (2 / 3) * mu * B * sqrt(2 * 9.81) * h^(3/2) +
    0.5 * (2 / 3) * mu * B * sqrt(2 * 9.81) * h^(3/2)
  Q <- uniroot(fcn, interval = c(1e-5, Q_max), check.conv = TRUE)$root
  v <- Q / ((h + w) * B)

  return(list(Q = Q, v = v))
}
