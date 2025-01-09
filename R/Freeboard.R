# -----------------------------------------------------------------------------
# Freeboard Calculation according to KOHS
# -----------------------------------------------------------------------------
#' @title Freeboard Calculation
#' @name freeboard
#' @description Calculates the required freeboard based on the KOHS (2013)
#' recommendations.
#' @usage freeboard(v, h, sigma_wz = 0, fw = TRUE, fv = FALSE, ft = 0, min = 0,
#'   max = Inf, fe_fixed = 0)
#' @param v Flow velocity [m/s].
#' @param h Flow depth [m].
#' @param sigma_wz Uncertainty in bed elevation (morphodynamics) [m].
#' @param fw Logical; considers freeboard due to uncertainty in water elevation.
#'   If `TRUE`, calculates according to KOHS; if `FALSE`, sets `fw = 0`.
#' @param fv Logical; considers freeboard due to waves. If `TRUE`, calculates
#'   according to KOHS; if `FALSE`, sets `fv = 0`.
#' @param ft Freeboard due to driftwood based on KOHS (2013) [m].
#' @param min Minimum allowable freeboard [m].
#' @param max Maximum allowable freeboard [m].
#' @param fe_fixed Fixed freeboard value to override calculations [m].
#' @return A numeric value of the calculated freeboard [m].
#' @references KOHS (2013). Freibord bei Hochwasserschutzprojekten und
#'   Gefahrenbeurteilungen - Empfehlungen der Kommission Hochwasserschutz KOHS.
#'   Wasser Energie Luft 105(1): 43-53.
#' @examples
#' freeboard(h = 1.36, sigma_wz = 0.3, fv = FALSE, ft = 0) # Channel example.
#' freeboard(v = 4.56, h = 1.36, sigma_wz = 0.3, fv = TRUE, ft = 0) # Dam.
#' freeboard(v = 4.56, h = 1.36, sigma_wz = 0.3, fv = TRUE, ft = 0.5) # Bridge.
#' @export


freeboard <- function(v, h, sigma_wz = 0, fw = TRUE, fv = FALSE,
                      ft = 0, min = 0, max = Inf, fe_fixed = 0) {
  if (fe_fixed == 0) {
    fw <- if (fw) sqrt(sigma_wz^2 + (0.06 + 0.06 * h)^2) else 0
    fv <- if (fv) v^2 / (2 * 9.81) else 0
    fe <- sqrt(fw^2 + fv^2 + ft^2)
    fe <- pmax(fe, min) # Constrain to minimum
    fe <- pmin(fe, max) # Constrain to maximum
  } else {
    fe <- fe_fixed
  }
  return(fe)
}
