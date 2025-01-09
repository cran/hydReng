#------------------------------------------------------------------------------
# Methods for class CScircle
#------------------------------------------------------------------------------

# Accessors
#------------------------------------------------------------------------------

setMethod("Di", "CScircle", function(object) object@Di)
setMethod("kSt", "CScircle", function(object) object@kSt)
setMethod("ks", "CScircle", function(object) object@ks)


# Validity
#------------------------------------------------------------------------------
setValidity("CScircle", function(object) {
  msg <- NULL
  valid <- TRUE
  if (Di(object) <= 0) {
    valid <- FALSE
    msg <- c(msg, "D must be a positive number.")
  }
  if (valid == TRUE) TRUE else msg
})



# Replacement
#------------------------------------------------------------------------------
setMethod("return_valid_object", "CScircle", function(object) {
  if (validObject(object) == TRUE) return(object)
}
)

setMethod("Di<-", "CScircle",
          function(object, value) {
            object@Di <- value
            return_valid_object(object)
          }
)

setMethod("kSt<-", "CScircle",
          function(object, value) {
            object@kSt <- value
            return_valid_object(object)
          }
)

setMethod("ks<-", "CScircle",
          function(object, value) {
            object@ks <- value
            return_valid_object(object)
          }
)


# geometry
#------------------------------------------------------------------------------

#wetted Area
#------------------------------------------------------------------------------
setMethod("wetted_area", "CScircle",
          function(object, h) {

            Di<-Di(object)
            if(Di>h){
              alpha<-2*acos(1-((Di-h)/(Di/2)))
              A_h<-((Di/2)^2*pi)-((Di/2)^2/2*(alpha-sin(alpha)))
            } else {A_h<-(Di/2)^2*pi}

            return(A_h)
          }
)


#wetted Perimeter
#------------------------------------------------------------------------------

setMethod("wetted_perimeter", "CScircle",
          function(object, h) {

            Di<-Di(object)
            if(Di>h){
              alpha<-2*acos(1-((Di-h)/(Di/2)))
              b<-(Di/2)*alpha
              P_h<-(Di*pi)-b
            } else {P_h<-Di*pi}

            return(P_h)
          }
)


# calculate velocity
#---------------------------------------------------
setMethod("flow_velocity", "CScircle",
          function(object, h, J, method, nu = 1.24e-6, ...) {

            if (method == "Strickler") {
              if (!is.numeric(kSt(object))) {
                warning("kSt is missing but mandatory for method 'Strickler'")
                return(NULL)
              }
              A <- wetted_area(object, h = h)
              P <- wetted_perimeter(object, h = h)
              v <- kSt(object) * sqrt(J) * (A / P)^(2 / 3)
            }

            else if (method == "Prandtl-Coolebrook-White") {
              if (!is.numeric(ks(object))) {
                warning("ks is missing but mandatory for method 'Prandtl-Coolebrook-White'")
                return(NULL)
              }
              A <- wetted_area(object, h = h)
              P <- wetted_perimeter(object, h = h)
              R <- A / P
              v <- (-2) * sqrt(8 * 9.81) * sqrt(R * J) * log10(
                ((ks(object) / 1000) / (14.8 * R)) +
                  ((2.51 * nu) / (4 * R * sqrt(8 * 9.81) * sqrt(R * J)))
              )
            }

            else {
              stop("Invalid method. Use 'Strickler' or 'Prandtl-Coolebrook-White'.")
            }

            return(v)
          }
)

# Froude number(Fr)
#---------------------------------------------------

setMethod("froude_number", "CScircle",
          function(object, v, h) {

            A <- wetted_area(object, h=h)
            Q <- v*A

            Fr<-Q/sqrt(9.81*Di(object)*h^4)
            return(Fr)

          }
)




# calculate flow depth
#---------------------------------------------------
setMethod("flow_depth", "CScircle",
          function(object, Q, J, method="Strickler", ret="all",plot=FALSE) {



            if(method=="Strickler"){

              if(!is.numeric(kSt(object))){
                warning("kSt is missing but mandatory for method 'Strickler'")
                return(NULL)
              }

              Qmax<-flow(object,h=Di(object)*0.9385,J=J,method=method)$Q
              if(Q<=Qmax){

                # function for searching root
                fcn <- function(h){
                  A <- wetted_area(object, h=h)
                  P <- wetted_perimeter(object, h=h)
                  v <- flow_velocity(object, h=h, J=J)
                  return(Q/A - v)
                }

                #search root
                if(is.numeric(try(uniroot(fcn, interval=c(10^-5, Di(object)*0.9385), check.conv = TRUE)$root,silent = TRUE))){
                  h <- uniroot(fcn, interval=c(10^-5, Di(object)*0.9385), check.conv = TRUE)$root
                  v<-flow_velocity(object, h=h, J=J,method=method)
                } else{
                  warning("h exceeds tube diameter")
                  h <- Di(object)
                  v<-flow_velocity(object, h=h, J=J,method=method)
                }

              } else {
                warning("h exceeds tube diameter")
                return(NULL)
              }
            }

            if(method=="Prandtl-Coolebrook-White"){           #Prandtl-Coolebrook-White according to SIA 190

              if(!is.numeric(ks(object))){
                warning("ks is missing but mandatory for method 'Prandtl-Coolebrook-White'")
                return(NULL)
              }

              Qmax<-flow(object,h=Di(object)*0.9385,J=J,method=method)$Q
              if(Q<=Qmax){

                # function for searching root
                fcn <- function(h){
                  A <- wetted_area(object, h=h)
                  P <- wetted_perimeter(object, h=h)
                  v <- flow_velocity(object, h=h, J=J,method=method)
                  return(Q/A - v)
                }

                #search root
                if(is.numeric(try(uniroot(fcn, interval=c(10^-5, Di(object)*0.9385), check.conv = TRUE)$root,silent = T))){
                  h <- uniroot(fcn, interval=c(10^-5, Di(object)*0.9385), check.conv = TRUE)$root
                  v<-flow_velocity(object, h=h, J=J,method=method)
                } else{
                  warning("h exceeds tube diameter")
                  h <- Di(object)
                  v<-flow_velocity(object, h=h, J=J,method=method)
                }

              } else{
                warning("h exceeds tube diameter")
                return(NULL)
              }
            }



            #plot

            if(plot==TRUE){
              try(dev.off(),silent=TRUE)

              if(h<=Di(object)){


                ylim<-if((h+ v^2/(2*9.81))<Di(object)){
                  c(0,Di(object))
                } else{c(0,h+ v^2/(2*9.81))}


                symbols(Di(object)/2,Di(object)/2,circles = Di(object)/2,inches=FALSE,xlim=c(0,Di(object)),ylim=ylim,asp=1, xlab="x [m]", ylab="z [m]")

                lines(c((Di(object)/2)-sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2),(Di(object)/2)+sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2)),c(h,h),col="blue")
                lines(c((Di(object)/2)-sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2),(Di(object)/2)+sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2)),rep(h+v^2/(2*9.81),2),col="red", lty=4)

                legend( "bottomleft",inset=c(0,1), xpd=TRUE,
                        legend=c("CS","water level","energy line","",
                                 paste("Q = ",round(Q,2)," m^3/s"),paste("J = ",J),paste("v = ",round(v,2)," m/s"),paste("F =",round(froude_number(object,v=v,h=h),2),", h=",round(h,2),"m")),
                        lty=c(1,1,2,NA,NA,NA,NA,NA),pch=c(NA,NA,NA,NA,NA,NA),col=c("black","blue","red","black",NA,NA,NA),
                        bty="n",cex=0.8,ncol=2
                )

              } else{
                warning("h exceeds tube diameter, no plot is drawn")
              }

            }

            #output
            if(ret=="all"){
              res <- list()
              res$h <- h
              res$v <- v
              res$Fr <-froude_number(object,v=v,h=h)
              res$A <- wetted_area(object, h=h)
              return(res)
            } else if (ret=="h"){
              return(h)
            } else if (ret=="v"){
              return(v)
            }

          }
)


# calculate flow (Q)
#---------------------------------------------------
setMethod("flow", "CScircle",
          function(object, h, J,method="Strickler", ret="all",plot=FALSE) {

            if(method=="Strickler"){

              if(!is.numeric(kSt(object))){
                warning("kSt is missing but mandatory for method 'Strickler'")
                return(NULL)
              }

              A <- wetted_area(object, h=h)
              v <- flow_velocity(object, h=h, J=J)
              Q <- v*A

            }

            if(method=="Prandtl-Coolebrook-White"){    #Prandtl-Coolebrook-White according to SIA 190

              if(!is.numeric(ks(object))){
                warning("ks is missing but mandatory for method 'Prandtl-Coolebrook-White'")
                return(NULL)
              }

              A <- wetted_area(object, h=h)
              U <- wetted_perimeter(object,h=h)
              R <- A/U
              v <- flow_velocity(object, h=h, J=J,method=method)
              Q <- v*A
            }



            #plot

            if(plot==TRUE){
              try(dev.off(),silent=TRUE)

              if(h<=Di(object)){


                ylim<-if((h+ v^2/(2*9.81))<Di(object)){
                  c(0,Di(object))
                } else{c(0,h+ v^2/(2*9.81))}


                symbols(Di(object)/2,Di(object)/2,circles = Di(object)/2,inches=FALSE,xlim=c(0,Di(object)),ylim=ylim,asp=1, xlab="x [m]", ylab="z [m]")

                lines(c((Di(object)/2)-sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2),(Di(object)/2)+sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2)),c(h,h),col="blue")
                lines(c((Di(object)/2)-sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2),(Di(object)/2)+sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2)),rep(h+v^2/(2*9.81),2),col="red", lty=4)

                legend( "bottomleft",inset=c(0,1), xpd=TRUE,
                        legend=c("CS","water level","energy line","",
                                 paste("Q = ",round(Q,2)," m^3/s"),paste("J = ",J),paste("v = ",round(v,2)," m/s"),paste("F =",round(froude_number(object,v=v,h=h),2),", h=",round(h,2),"m")),
                        lty=c(1,1,2,NA,NA,NA,NA,NA),pch=c(NA,NA,NA,NA,NA,NA),col=c("black","blue","red","black",NA,NA,NA),
                        bty="n",cex=0.8,ncol=2
                )

              } else{
                warning("h exceeds tube diameter, no plot is drawn")
              }

            }



            #output
            if(ret=="all"){
              res <- list()
              res$Q <- Q
              res$v <- v
              res$A <- A
              res$Fr<-froude_number(object,v=v,h=h)
              return(res)
            } else if (ret=="Q"){
              return(Q)
            } else if (ret=="v"){
              return(v)
            }

          }
)

# flow_max
#---------------------------------------------------


setMethod("flow_max", "CScircle",
          function(object, J,method="Strickler", ret="all",plot=FALSE) {


            if(method=="Strickler"){

              if(!is.numeric(kSt(object))){
                warning("kSt is missing but mandatory for method 'Strickler'")
                return(NULL)
              }

              #flow for maximal level
              v.max<-c()
              for(i in seq(1,(Di(object)*100),1)){
                v.max[i]<-flow(object,i/100,J=J,ret="Q")
              }
              Qmax <- max(v.max)
              hmax<-which.max(v.max)/100
              v <- flow_velocity(object, h=hmax, J=J)

            }


            if(method=="Prandtl-Coolebrook-White"){    #Prandtl-Coolebrook-White according to SIA 190

              if(!is.numeric(ks(object))){
                warning("ks is missing but mandatory for method 'Prandtl-Coolebrook-White'")
                return(NULL)
              }

              #flow for maximal level
              v.max<-c()
              for(i in seq(1,(Di(object)*100),1)){
                v.max[i]<-flow(object,i/100,J=J,method=method, ret="Q")
              }
              Qmax <- max(v.max)
              hmax<-which.max(v.max)/100
              v <- flow_velocity(object, h=hmax, J=J,method=method)

            }


            #plot
            h<-hmax
            Q<-Qmax

            if(plot==TRUE){
              try(dev.off(),silent=TRUE)


              ylim<-if((h+ v^2/(2*9.81))<Di(object)){
                c(0,Di(object))
              } else{c(0,h+ v^2/(2*9.81))}


              symbols(Di(object)/2,Di(object)/2,circles = Di(object)/2,inches=FALSE,xlim=c(0,Di(object)),ylim=ylim,asp=1, xlab="x [m]", ylab="z [m]")

              lines(c((Di(object)/2)-sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2),(Di(object)/2)+sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2)),c(h,h),col="blue")
              lines(c((Di(object)/2)-sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2),(Di(object)/2)+sqrt((Di(object)/2)^2-(abs(Di(object)/2-h))^2)),rep(h+v^2/(2*9.81),2),col="red", lty=4)

              legend( "bottomleft",inset=c(0,1), xpd=TRUE,
                      legend=c("CS","water level","energy line","",
                               paste("Q = ",round(Q,2)," m^3/s"),paste("J = ",J),paste("v = ",round(v,2)," m/s"),paste("F =",round(froude_number(object,v=v,h=h),2),", h=",round(h,2),"m")),
                      lty=c(1,1,2,NA,NA,NA,NA,NA),pch=c(NA,NA,NA,NA,NA,NA),col=c("black","blue","red","black",NA,NA,NA),
                      bty="n",cex=0.8,ncol=2
              )

            }

            #output
            if(ret=="all"){
              res <- list()
              res$Qmax <- Qmax
              res$hmax <- hmax
              res$v <- flow_velocity(object, h=hmax, J=J,method=method)
              res$A <- wetted_area(object, h=hmax)
              res$Fr<-froude_number(object,v=v,h=hmax)
              return(res)
            } else if (ret=="Qmax"){
              return(Qmax)
            } else if (ret=="hmax"){
              return(hmax)
            }else if (ret=="v"){
              return(flow_velocity(object, h=hmax, J=J,method=method))
            }

          }
)


# Partial filling diagram
#---------------------------------------------------
#' @title Partial Filling Flow Diagram
#' @name par_fill
#' @description Function to generate a plot of partial-filling diagram of
#'  circular pipe with discharge and flow velocity
#' @aliases par_fill,CScircle-method
#' @usage par_fill(object,J,method="Strickler")
#' @param object A CScircle object.
#' @param J Bottom slope [-].
#' @param method Method to calculate the roughness. Allowed are "Strickler"
#'  (equal roughness) and "Prandtl-Coolebrook-White".
#' @return Plots of a partial filling diagram of a circular pipe with discharge and flow velocity
#' @examples
#' csC <- CScircle(Di = 0.7, ks = 1.5, kSt = 75)
#' par_fill(csC,J=0.04)
#' @export
#'
#'


setMethod("par_fill", "CScircle",
          function(object,J,method="Strickler"){

            if(method=="Strickler"){

              if(!is.numeric(kSt(object))){
                warning("kSt is missing but mandatory for method 'Strickler'")
                return(NULL)
              }

              v.vel<-c()
              v.Q<-c()
              for(i in 1: (Di(object)*100)){
                v.vel[i]<-flow(object,h=i/100,J=J)$v
                v.Q[i]<-flow(object,h=i/100,J=J)$Q

              }

            }

            if(method=="Prandtl-Coolebrook-White"){

              if(!is.numeric(ks(object))){
                warning("ks is missing but mandatory for method 'Prandtl-Coolebrook-White'")
                return(NULL)
              }

              v.vel<-c()
              v.Q<-c()
              for(i in 1: (Di(object)*100)){
                v.vel[i]<-flow(object,h=i/100,J=J,method=method)$v
                v.Q[i]<-flow(object,h=i/100,J=J,method=method)$Q

              }

            }

            plot(1,1,xlim=c(0,max(max(v.Q),max(v.vel))),ylim=c(0,(Di(object)*100)),type="n",yaxt="n",ylab="Filling [%]",xlab="Q [m3/s]; v [m/s]",main="Partial filling - flow diagram")
            points((v.Q),1:(Di(object)*100),t="l")
            points((v.vel),1:(Di(object)*100),t="l",lty=2)
            #axis(1,at=seq(0,max(max(v.vel),max(v.Q)),0.5),labels=seq(0,max(max(v.vel),max(v.Q)),0.5))
            axis(2,at=seq(0,Di(object)*100,length.out = 5),labels=c(0,25,50,75,100))
            legend("bottomright",legend=c("Discharge","Velocity"),lty=c(1,2))
          }
)
