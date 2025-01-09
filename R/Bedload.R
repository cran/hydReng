#------------------------------------------------------------------------------
## Bedload transport
#------------------------------------------------------------------------------

#' @title Bedload Transport Capacity (Smart and Jaeggi)
#' @description
#' Calculates the bedload transport capacity based on the formula by Smart and
#' Jaeggi (1983). This formula is recommended for slopes between 0.005 and 0.2.
#'
#' @usage
#' bedload_SJ(d30, dm, d90, J, Rs, um, B, t_crit = 0.05, rho_s = 2650,
#' s_value = 2.65)
#'
#' @param d30 Grain size distribution parameter [m].
#' @param dm Median grain size [m].
#' @param d90 Grain size distribution parameter [m].
#' @param J Bottom slope [-].
#' @param Rs Hydraulic radius [m].
#' @param um Mean flow velocity [m/s].
#' @param B Bottom width [m].
#' @param t_crit Critical shear stress [-] (default: 0.05).
#' @param rho_s Density of bedload material [kg/m3] (default: 2650).
#' @param s_value Relative solid density [-] (default: 2.65).
#'
#' @return
#' bedload_SJ returns the bedload transport rate [kg/s]
#'
#'
#' @references
#' Smart, G. M., & Jäggi, M. N. R. (1983). Sediment transport in steilen
#' Gerinnen. Mitteilungen der Versuchsanstalt für Wasserbau, Hydrologie und
#' Glaziologie der ETH Zürich, 64, Zürich.
#'
#' @examples
#' d30 <- 0.05
#' dm <- 0.1
#' d90 <- 0.2
#' J <- 0.03
#' Rs <- 1
#' um <- 2
#' B <- 3
#'
#' bedload_SJ(d30 = 0.05, dm = 0.10, d90 = 0.2, J = 0.03, Rs = 1, um = 2, B = 5)
#'
#' @export



bedload_SJ <- function(d30, dm, d90, J, Rs, um, B, t_crit = 0.05, rho_s = 2650,
                  s_value = 2.65) {

  Gb <- max(0, 4 * rho_s / (s_value - 1) * (d90 / d30)^0.2 * J^1.6 * Rs * um *
              B * (1 - t_crit * (s_value - 1) * dm / (Rs * J)))

  return(Gb)
}


#' @title Bedload Transport Capacity (Meyer-Peter Müller)
#' @description
#' Calculates the bedload transport capacity using the formula by
#' Meyer-Peter Müller. The formula is valid for bed slopes less than 0.005.
#'
#' @usage
#' bedload_MPM(dm, J, Rs, B, f_kSt = 0.85, t_crit = 0.047, rho_s = 2650, s = 2.65)
#'
#' @param dm Median grain size [m].
#' @param J Bottom slope [-].
#' @param Rs Hydraulic radius [m].
#' @param B Bottom width [m].
#' @param f_kSt Friction factor = (k_StS / k_Str)^(3/2) (default: 0.85).
#' @param t_crit Critical shear stress [-] (default: 0.047).
#' @param rho_s Density of bedload material [kg/m3] (default: 2650).
#' @param s Relative solid density [-] (default: 2.65).
#'
#' @return
#' Returns the bedload transport rate [kg/s].
#'
#' @references
#' Bezzola, G.R. (2012). Vorlesungsmanuskript Flussbau. ETH Zürich, Versuchsanstalt für Wasserbau, Hydrologie und Glaziologie VAW.
#'
#' @examples
#' bedload_MPM(dm = 0.1, J = 0.01, Rs = 1.5, B = 20)
#' bedload_MPM(dm = 0.1, J = 0.01, Rs = 1.5, B = 20, t_crit = 0.06)
#'
#' @export





bedload_MPM <- function(dm, J, Rs, B, f_kSt=0.85,t_crit = 0.047, rho_s=2650,s=2.65){

    if(f_kSt*Rs*J - t_crit*(s-1)*dm > 0){
      Gb<- ((8*sqrt(9.81)*rho_s*B)/(s-1)) * (f_kSt*Rs*J - t_crit*(s-1)*dm)^(3/2)
    } else {
      Gb <- 0
    }


  return(Gb)

}
