# ------------------------------------------------------------------------------
## vertical scour
# ------------------------------------------------------------------------------

#' Scour depth and length (vertical jet)
#'
#' Calculate scour depth and position (length) formed by vertical jet
#'
#' @usage
#' scour_vert(
#'   Q, B, h, h_u, d90, d95, ful_ov, method = "Kolatus", bedload = FALSE
#' )
#'
#' @param Q
#' flow [m3/s]
#'
#' @param B
#' width of the overfall section [m]
#'
#' @param h
#' difference between upstream and downstream water table [m]
#'
#' @param h_u
#' downstream water table [m]
#'
#' @param d90
#' d90 of grain size distribution [mm]
#'
#' @param d95
#' d95 of grain size distribution [mm]
#'
#' @param ful_ov
#' defines if the overfall is complete or incomplete. TRUE and FALSE are valid values
#' [logical]
#'
#' @param method
#' method to calculate scour depth. valid values are "Kolatus" or "Tschopp".
#' Independent of the chosen method, the position of the scour is calculated
#' according to "Kolatus" and the scour depth considering bedload is calculated
#' according to "Tschopp".
#'
#' @param bedload
#' Consider bedload transportation. If bedload = TRUE, the scour depth "TG" and
#' "SG" are calculated additionally which consider bedload transport.
#'
#' @return
#' \item{T0}{water table at maximal scour depth [m]}
#' \item{S}{maximal scour depth [m]}
#' \item{TG}{water table at maximal scour depth considering bedload transport
#'   [m]}
#' \item{SG}{maximal scour depth considering bedload transport [m]}
#' \item{l1}{horizontal distance of maximal scour depth to overfall crest [m]}
#' \item{l2}{total horizontal distance of scour depth from overfall crest [m]}
#'
#' @references
#' Bezzola (2012). Vorlesungsmanuskript Flussbau. ETH Zürich.
#'
#' @examples
#' ## calculate scour depth according to Kolatus returning all results
#' scour_vert(
#'   Q = 4, B = 1, h = 3, h_u = 1.76, d90 = 150, d95 = 200, ful_ov = TRUE
#' )
#'
#' ## calculate scour depth according to Tschopp
#' scour_vert(
#'   Q = 4, B = 1, h = 3, h_u = 1.76, d90 = 150, d95 = 200,
#'   method = "Tschopp", ful_ov = TRUE
#' )$S
#'
#' ## calculate scour depth according to Tschopp considering bedload transport
#' scour_vert(
#'   Q = 4, B = 1, h = 3, h_u = 1.76, d90 = 150, d95 = 200,
#'   method = "Tschopp", bedload = TRUE, ful_ov = TRUE
#' )$SG
#'
#' @export


scour_vert <- function(Q, B, h, h_u, d90, d95, ful_ov,
                       method = "Kolatus", bedload = FALSE) {

  q <- Q / B

  if (ful_ov == TRUE) { # full Overfall

    if (
      ((q^(1 / 2) * h^(1 / 4)) / (d90 / 1000) >= 5 &
       (q^(1 / 2) * h^(1 / 4)) / (d90 / 1000) <= 25) == FALSE
    ) { # Gueltigkeitsbereich
      warning("validity criterion is not fulfilled")
    }

    if (method == "Kolatus") {
      T0 <- 0.78 * (h^0.35 * q^0.7) / ((d90 / 1000)^0.4)
      S <- T0 - h_u
      l1 <- 3.9 * (
        (h^0.27 * q^0.54) /
          (9.81^0.27 * (d95 / 1000)^0.08)
      )
      l2 <- 2.7 * (
        (h^0.45 * q^0.9) /
          (9.81^0.45 * (d95 / 1000)^0.8)
      )
    }

    if (method == "Tschopp") {
      T0 <- 2.76 * q^(1 / 2) * h^(1 / 4) - 7.22 * (d90 / 1000)
      S <- T0 - h_u

      if (bedload == TRUE) {
        TG <- 1.9 * q^0.5 * h^0.25 - 4.17 * (d90 / 1000)
        SG <- TG - h_u

        if (
          (5 <= TG / (d90 / 1000) &
           60 >= TG / (d90 / 1000) &
           1.1 <= TG / (q^0.5 * h^0.25) &
           1.8 >= TG / (q^0.5 * h^0.25)) == FALSE
        ) {
          warning(
            "validity criterion for beadload transport is not fulfilled"
          )
        }
      }

      l1 <- 3.9 * (
        (h^0.27 * q^0.54) /
          (9.81^0.27 * (d95 / 1000)^0.08)
      )
      l2 <- 2.7 * (
        (h^0.45 * q^0.9) /
          (9.81^0.45 * (d95 / 1000)^0.8)
      )
    }

  } else { # unvollkommener Überfall

    if (method == "Kolatus") {
      T0 <- (0.78 * (h^0.35 * q^0.7) / ((d90 / 1000)^0.4)) * 0.84
      S <- T0 - h_u
      l1 <- 3.9 * (
        (h^0.27 * q^0.54) /
          (9.81^0.27 * (d95 / 1000)^0.08)
      ) * 1.13
      l2 <- 2.7 * (
        (h^0.45 * q^0.9) /
          (9.81^0.45 * (d95 / 1000)^0.8)
      ) * 1.6
    }

    if (method == "Tschopp") {
      T0 <- (2.76 * q^(1 / 2) * h^(1 / 4) - 7.22 * (d90 / 1000)) * 0.84
      S <- T0 - h_u

      if (bedload == TRUE) {
        TG <- 1.9 * q^0.5 * h^0.25 - 4.17 * (d90 / 1000)
        SG <- TG - h_u

        if (
          (5 <= TG / (d90 / 1000) &
           60 >= TG / (d90 / 1000) &
           1.1 <= TG / (q^0.5 * h^0.25) &
           1.8 >= TG / (q^0.5 * h^0.25)) == FALSE
        ) {
          warning(
            "validity criterion for beadload transport is not fulfilled"
          )
        }
      }

      l1 <- 3.9 * (
        (h^0.27 * q^0.54) /
          (9.81^0.27 * (d95 / 1000)^0.08)
      ) * 1.13
      l2 <- 2.7 * (
        (h^0.45 * q^0.9) /
          (9.81^0.45 * (d95 / 1000)^0.8)
      ) * 1.6
    }
  }

  res <- list()

  if (bedload == FALSE) {
    res$T0 <- T0
    res$S <- S
    res$TG <- NA
    res$SG <- NA
    res$l1 <- l1
    res$l2 <- l2
  }

  if (bedload == TRUE) {
    res$T0 <- T0
    res$S <- S
    res$TG <- TG
    res$SG <- SG
    res$l1 <- l1
    res$l2 <- l2
  }

  return(res)
}



# horizontal scour
#---------------------------------------------------

#' Scour depth and length (horizontal jet)
#'
#' Calculate scour depth and position (length) formed by horizontal jet
#'
#' @usage scour_horz(Q,B,h,h_u,d90,W=NULL,a=NULL,mu=NULL,l0=NULL,
#'method="Eggenberger")

#'
#' @param Q flow [m3/s]
#' @param B channel width [m]
#' @param h difference between upstream and downstream water table [m]
#' @param h_u downstream water table [m]
#' @param d90 d90 of grain size distribution [mm]
#' @param W shape value. For free waved jet w = 15.4 and for covered waved jet
#'   w = 10.35. If the parameters a and mu are known, W is calculated in the
#'   function and must not be specified
#' @param a gate opeing height [m]
#' @param mu contraction value. Values between 0.6 (orthogonal gates) and 0.8
#'   (inclined gates) are recommended [-]
#' @param l0 length of fix sole protection downstream the gate [m] (for method Shalash)
#' @param method method to calculate scour depth. valid values are "Eggenberger"
#'   or "Shalash". Independet of the chosen method, the position of the scour is
#'   calculated according to "Eggenberger".
#'
#'
#' @return
#' \item{T0}{water table at maximal scour depth [m]}
#' \item{S}{maximal scour depth [m]}
#' \item{l1}{horizontal distance of maximal scour depth to overfall crest [m]}
#' \item{l2}{total horizontal distance of scour depth from overfall crest [m]}
#'
#' @references
#' Bezzola (2012). Vorlesungsmanuskript Flussbau. ETH Zürich.
#'
#'
#' @examples
#' ## calculate scour depth accordint to Eggenberger returing all results
#' scour_horz(Q = 4, B = 1, h = 5, h_u = 1.76, d90 = 150, a = 1, mu = 0.6)
#'
#' ## calculate scour depth accordint to Eggenberger with given W value
#' scour_horz(Q = 4, B = 1, h = 5, h_u = 1.76, d90 = 150, W = 15.4)
#'
#' ## calculate scour depth accordint to Shalash with l0=5
#' scour_horz(
#'   Q = 4, B = 1, h = 5, h_u = 1.76, d90 = 150,
#'   method = "Shalash", l0 = 5, a = 1
#' )
#'
#' @export




scour_horz<-function(Q,B,h,h_u,d90,W=NULL,a=NULL,mu=NULL,l0=NULL,method="Eggenberger"){

  q<-Q/B

  if(method=="Eggenberger"){

    if(is.null(W)){

      if(any(is.null(a),is.null(mu))){
        warning("if W is NULL, a and mu must not be equal NULL")
      }


    if(h_u<=-(mu*a/2)+sqrt(((mu^2*a^2)/4)+(2*q^2/(9.81*mu*a)))){  #freier Durchfluss
      W<-15.4
    } else {    #bedeckter Durchfluss
      W<-10.35
    }
    }


    T0<-W*(h^0.5*q^0.6)/d90^0.4
    S<-T0-h_u
  }

if(method=="Shalash"){
  T0<-9.65*((h^0.5*q^0.6)/d90^0.4)*(1.5*h/l0)^0.6
  S<-T0-h_u

}


  l1<-4.9*S
  l2<-9.9*S

  res<-list()
  res$T0<-T0
  res$S<-S
  res$l1<-l1
  res$l2<-l2

  return(res)

}


# scour in curve
#---------------------------------------------------

#' Scour depth in a curve
#'
#' Calculate scour depth formed in a curve.
#'
#' @usage scour_curve(A,v,Fr,h,J,r,rm,d84,d16,dm=NULL,psi=NULL,method="Peter")
#'
#' @param A wetted area upstream the curve [m2]
#' @param v flow velocity upstream the curve [m/s]
#' @param Fr Froude number upstream the curve [-]
#' @param h flow depth in the middle of the river upstream the curve [m]
#' @param J bottom slope [-]
#' @param r curve radius from center to the outer bank bottom [m]
#' @param rm curve radius from center to the middle of the river [m]
#' @param d16 d16 of grain size distribution [mm]
#' @param dm d50 of grain size distribution [mm]
#' @param d84 d84 of grain size distribution [mm]
#' @param psi inner friction angle[°]. Values between 20° and 25° are
#'   recommended for flat rivers (J~0.0003). For steeper rivers (0.0035 < J <
#'   0.007), values between 35° and 40° are recomended (Bezzola 2012).
#' @param method method to calculate scour depth. valid values are "Peter",
#'   "Bridge" or "Kikkawa"
#'
#' @return
#' \item{T0}{water table at maximal scour depth [m]}
#' \item{S}{difference between bed elevation at the middle of the river upstream
#'   the curve and the maximal scour depth  [m]}
#' \item{val}{if val = T, the criterion for method "Peter" is fullfilled}
#'
#' @references
#' Bezzola (2012). Vorlesungsmanuskript Flussbau. ETH Zürich.
#'
#'
#' @examples
#' ## calculate scour depth accordint to Peter
#'
#' # Define parameter
#' A <- 135.5
#' Fr <- 0.52
#' h <- 3.31
#' J <- 0.0022
#' r <- 500
#' rm <- 530
#' d16 <- 50
#' d84 <- 200
#'
#' scour_curve(
#'   A = A, Fr = Fr, h = h, J = J, rm = rm, r = r,
#'   d16 = d16, d84 = d84
#' )
#'
#' ## calculate scour depth accordint to Bridge
#'
#' # Define parameter
#' h <- 3.31
#' r <- 500
#' rm <- 530
#' psi <- 20
#'
#' scour_curve(h = h, rm = rm, r = r, method = "Bridge", psi = psi)
#'
#' ## calculate scour depth according to Kikkawa
#'
#' # Define parameter
#' v <- 2.7
#' h <- 3.31
#' J <- 0.0022
#' r <- 500
#' rm <- 470
#' dm <- 80
#'
#' scour_curve(
#'   v = v, h = h, J = J, rm = rm, r = r,
#'   dm = dm, method = "Kikkawa"
#' )
#'
#' @export
scour_curve <- function(A, v, Fr, h, J, r, rm, d16,
                        dm = NULL, d84, psi = NULL,
                        method = "Peter") {
  # function body
}


scour_curve<-function(A,v,Fr,h,J,r,rm,d84,d16,dm=NULL,psi=NULL,method="Peter"){

  if(method=="Peter"){
    bm<-A/h
    sigma<-sqrt(d84/d16)


    K<-2.95*(rm/bm)-0.7*sigma-29.3*(h/bm)+2.7*Fr+3.4

    if((rm/bm)>=2 & (rm/bm)<=6){
      val<-TRUE
    }else{
      val<-FALSE
    }
    T0<-h*(r/rm)^K
    S<-T0-h


  }


  if(method=="Bridge"){

    if(!is.numeric(psi)){
      warning("psi must be numeric")
      return(NULL)
    }


    K<-11*tan(psi*(pi/180))
    T0<-h*(r/rm)^K
    S<-T0-h

    val<-NA

  }

  if(method=="Kikkawa"){

    if(!is.numeric(dm)){
      warning("dm must be numeric")
      return(NULL)
    }

    nu <- sqrt(shear_str(h0=h,J=J,dm=dm)$tau_) * ( 1.9 * v/shear_str(h0=h,J=J)$U - 3 )
    T0 <- h * exp(0.5*nu*(r^2/rm^2-1))
    S<-T0-h

    val<-NA

  }

  return(list(T0=T0,S=S,val=val))
}







# scour at groyne
#---------------------------------------------------

#' Scour depth formed by a groyne
#'
#' Calculate scour depth formed by a groyne
#'
#' @usage scour_groyne(v,Fr,B,h,J,L,d16,dm,d84,Ks,delta,Kb=NULL,l=NULL, fs=0,
#' method="Froehlich",bedload=FALSE)

#' @param v flow velocity upstream the groyne [m/s]
#' @param Fr Froude number upstream the groyne [-]
#' @param B sole width [m]
#' @param h flow depth upstream the groyne [m]
#' @param J bottom slope [-]
#' @param L length of the groyne (perpendicular to the river) [m]
#' @param d16 d16 of grain size distribution [mm]
#' @param dm d50 of grain size distribution [mm]
#' @param d84 d84 of grain size distribution [mm]
#' @param Ks shape value according to "Froehlich". values between 0.55 and 1
#'   are recomended.[-]
#' @param delta horizontal angle of the groyne in respect to the river [°]
#' @param Kb shape value according to "Hoffmanns" [-]
#' @param l length of the groyne (parallel to the river)[m]
#' @param fs safety factor [-]
#' @param method method to calculate scour depth. valid values are "Froehlich"
#' @param bedload Consider bedload transportation if bedload =TRUE
#'
#' @return
#' \item{T0}{water table at maximal scour depth [m]}
#' \item{S}{difference between bed elevation at the middle of the river and
#'   the maximal scour depth  [m]}
#'
#' @references
#' Bezzola (2012). Vorlesungsmanuskript Flussbau. ETH Zürich.
#'
#'
#' @examples
#' ## calculate scour depth accordint to Froehlich without bedload
#'
#' v <- 2.7
#' Fr <- 0.52
#' h <- 3.31
#' J <- 0.0022
#' L <- 5
#' d16 <- 50
#' dm <- 80
#' d84 <- 200
#' Ks <- 0.82
#' delta <- 60
#'
#' scour_groyne(
#'   v = v, Fr = Fr, h = h, J = J, L = L,
#'   d16 = d16, dm = dm, d84 = d84,
#'   Ks = Ks, delta = delta
#' )
#'
#' ## calculate scour depth accordint to Froehlich with bedload
#'
#' v <- 2.7
#' Fr <- 0.52
#' h <- 3.31
#' J <- 0.0022
#' L <- 5
#' Ks <- 0.82
#' delta <- 60
#'
#' scour_groyne(
#'   v = v, Fr = Fr, h = h, J = J, L = L,
#'   d16 = d16, dm = dm, d84 = d84,
#'   Ks = Ks, delta = delta, bedload = TRUE
#' )
#'
#' @export



scour_groyne <- function(v,Fr,B,h,J,L,d16,dm,d84,Ks,delta,Kb=NULL,l=NULL, fs=0,
                         method="Froehlich",bedload= FALSE){
  if(method=="Froehlich"){

    K_delta<-(delta/90)^0.1

    if(bedload==FALSE){
      sigma<-sqrt(d84/d16)
      S<-(0.78*Ks*K_delta*(L/h)^0.63*Fr*(h/(dm/1000))^0.43*sigma^-1.87*h)+fs
    } else{
      S<-(2.27*K_delta*Ks*(L/h)^0.43*Fr^0.61*h)+fs

    }
    T0<-S+h


  }


  # if(method=="Hoffmanns"){
  #
  #   if(bedload==FALSE){
  #     warning("bedload must be TRUE")
  #     return(NULL)
  #   }
  #
  #   if(!is.numeric(Kb)){
  #     warning("Kb must be numeric")
  #     return(NULL)
  #   }
  #
  #
  #   m<-(B-L)/L
  #
  #   S<-h*((l-m)^-(2/3)-1)+Kb*L*tanh(h/L)
  #   T0<-S+h
  #
  # }

  return(list(T0=T0,S=S))

}










