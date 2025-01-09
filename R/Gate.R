#------------------------------------------------------------------------------
# Functions for calculating flow under gate
#------------------------------------------------------------------------------
#--------------------------------------------------------------
# psi_gate --> Bollrich 5. Aufl. (8.25) und (8.25a)
#--------------------------------------------------------------

psi_gate <- function(a, h0, alpha = 90) {
  if (alpha == 90) { # Vertical gate
    psi <- 1 / (1 + 0.64 * sqrt(1 - (a / h0)^2))
  } else if (alpha < 90 & alpha >= 30) { # Inclined gate
    psi <- 1.3 - 0.8 * sqrt(1 - ((alpha - 205) / 220)^2)
  } else {
    stop("Alpha must be between 30 and 90 degrees.")
  }
  return(psi)
}

#--------------------------------------------------------------
# mu_gate --> Bollrich 5. Aufl. (8.23)
#--------------------------------------------------------------

mu_gate <- function(psi, a, h0) {
  mu <- psi / sqrt(1 + (psi * a) / h0)
  return(mu)
}

#--------------------------------------------------------------
# flow_gate --> Bollrich 5. Aufl. (8.23)
#--------------------------------------------------------------

#' @title Discharge At Underflow Gate
#' @description Calculates the discharge through a gate under free or submerged
#'   conditions.
#' @usage flow_gate(a, h0, B, alpha, h2 = NULL, ret = "Q")
#' @param a Gate opening height [m].
#' @param h0 Upstream water depth [m].
#' @param B Gate width [m].
#' @param alpha Gate angle from horizontal [degrees].
#' @param h2 Optional. Downstream water depth [m]. Default is NULL (free flow).
#' @param ret Specifies the return value. "Q" for discharge only or "all" for
#'   all intermediate results.
#' @return A list containing the following hydraulic variables:
#' \describe{
#'   \item{Q}{Flow [m3/s].}
#'   \item{psi}{Contraction coefficient [-].}
#'   \item{mu}{Discharge coefficient [-].}
#'   \item{v}{Flow velocity [m/s].}
#'   \item{chi}{Coefficient for submerged flow [-].}
#' }
#' @examples
#' flow_gate(a = 0.5, h0 = 1.0, B = 2.0, alpha = 90)
#' flow_gate(a = 0.5, h0 = 1.0, B = 2.0, alpha = 90, h2 = 0.8)
#' flow_gate(a = 0.5, h0 = 1.0, B = 2.0, alpha = 90, h2 = 0.8, ret = "all")
#' @export
flow_gate <- function(a, h0, B, alpha, h2 = NULL, ret = "Q") {
  psi <- psi_gate(a, h0, alpha)
  mu <- mu_gate(psi, a, h0)
  v <- sqrt(2 * 9.81 * h0)

  if (is.null(h2)) { # Free flow
    Q <- mu * B * a * sqrt(2 * 9.81 * h0)
    chi <- NA
  } else { # Submerged flow
    expr <- 1 - (2 * psi * a / h0) * (1 - psi * a / h2)
    chi <- ifelse((expr^2 + (h2 / h0)^2 - 1) > 0,
                  sqrt((1 + psi * a / h0) * (expr - sqrt(expr^2 + (h2 / h0)^2 - 1))),
                  1)
    Q <- chi * mu * B * a * sqrt(2 * 9.81 * h0)
  }

  if (ret == "Q") {
    return(Q)
  } else if (ret == "all") {
    return(list(psi = psi, mu = mu, Q = Q, v = v, chi = chi))
  } else {
    stop("Invalid return type. Use 'Q' or 'all'.")
  }
}

#--------------------------------------------------------------
# flow_depth_gate --> Bollrich 5. Aufl. (8.23)
#--------------------------------------------------------------

#' @title Water Depth Upstream Of Gate
#' @description Calculates the upstream water depth for a gate based on given
#'   discharge and gate parameters.
#' @usage flow_depth_gate(a, Q, B, alpha, h2 = NULL, ret = "h0")
#' @param a Gate opening height [m].
#' @param Q Discharge [m3/s].
#' @param B Gate width [m].
#' @param alpha Gate angle from horizontal [degrees].
#' @param h2 Optional. Downstream water depth [m]. Default is NULL (free flow).
#' @param ret Specifies the return value. "h0" for depth only or "all" for
#'   all intermediate results.
#' @return A list containing the following hydraulic variables:
#' \describe{
#'   \item{h0}{Upstream water depth [m].}
#'   \item{psi}{Contraction coefficient [-].}
#'   \item{mu}{Discharge coefficient [-].}
#'   \item{v}{Flow velocity [m/s].}
#' }
#' @examples
#' flow_depth_gate(a = 0.5, Q = 2.5, B = 2.0, alpha = 90)
#' flow_depth_gate(a = 0.5, Q = 2.5, B = 2.0, alpha = 90, h2 = 0.8)
#' flow_depth_gate(a = 0.5, Q = 2.5, B = 2.0, alpha = 90, h2 = 0.8, ret = "all")
#' @export
flow_depth_gate <- function(a, Q, B, alpha, h2 = NULL, ret = "h0") {
  fcn <- function(h0) Q - flow_gate(a = a, h0 = h0, B = B, alpha = alpha, h2 = h2, ret = "Q")
  low_interval <- if (is.null(h2)) 1 else h2 + 0.001
  h0 <- uniroot(fcn, interval = c(low_interval, 10000), check.conv = TRUE)$root

  psi <- psi_gate(a, h0, alpha)
  mu <- mu_gate(psi, a, h0)
  v <- sqrt(2 * 9.81 * h0)

  if (ret == "h0") {
    return(h0)
  } else if (ret == "all") {
    return(list(psi = psi, mu = mu, h0 = h0, v = v))
  } else {
    stop("Invalid return type. Use 'h0' or 'all'.")
  }
}
