#------------------------------------------------------------------------------
# Methods for class CSarbitrary
#------------------------------------------------------------------------------

# Accessors
#------------------------------------------------------------------------------

setMethod("x", "CSarbitrary", function(object) object@x)
setMethod("z", "CSarbitrary", function(object) object@z)
setMethod("xb_l", "CSarbitrary", function(object) object@xb_l)
setMethod("xb_r", "CSarbitrary", function(object) object@xb_r)
setMethod("kSt_B", "CSarbitrary", function(object) object@kSt_B)
setMethod("kSt_l", "CSarbitrary", function(object) object@kSt_l)
setMethod("kSt_r", "CSarbitrary", function(object) object@kSt_r)

# Validity check
#------------------------------------------------------------------------------

setValidity("CSarbitrary", function(object) {
  msg <- NULL
  valid <- TRUE

  # Check if lengths of x and z match
  if (length(x(object)) != length(z(object))) {
    valid <- FALSE
    msg <- c(msg, "x and z have not the same length")
  }

  # Check xb_l if it's not NULL and if it matches any x
  if (!is.null(xb_l(object))) {
    if (!(xb_l(object) %in% x(object))) {
      valid <- FALSE
      msg <- c(msg, "xb_l is not equal to any x")
    }
  }

  # Check xb_r if it's not NULL and if it matches any x
  if (!is.null(xb_r(object))) {
    if (!(xb_r(object) %in% x(object))) {
      valid <- FALSE
      msg <- c(msg, "xb_r is not equal to any x")
    }
  }

  # Check if x is unique
  if (any(duplicated(x(object)))) {
    valid <- FALSE
    msg <- c(msg, "x values must be unique")
  }

  # Return validation result
  if (valid) TRUE else msg
})

# Replacement methods
#------------------------------------------------------------------------------
setMethod("return_valid_object", "CSarbitrary", function(object) {
  if (validObject(object)) return(object)
})

setMethod("x<-", "CSarbitrary", function(object, value) {
  object@x <- value
  return_valid_object(object)
})

setMethod("z<-", "CSarbitrary", function(object, value) {
  object@z <- value
  return_valid_object(object)
})

setMethod("xb_r<-", "CSarbitrary", function(object, value) {
  object@xb_r <- value
  return_valid_object(object)
})

setMethod("xb_l<-", "CSarbitrary", function(object, value) {
  object@xb_l <- value
  return_valid_object(object)
})

setMethod("kSt_B<-", "CSarbitrary", function(object, value) {
  object@kSt_B <- value
  return_valid_object(object)
})

setMethod("kSt_l<-", "CSarbitrary", function(object, value) {
  object@kSt_l <- value
  return_valid_object(object)
})

setMethod("kSt_r<-", "CSarbitrary", function(object, value) {
  object@kSt_r <- value
  return_valid_object(object)
})

# Wetted Area
#------------------------------------------------------------------------------
#' @title Wetted Area
#' @name wetted_area
#' @aliases wetted_area wetted_area,CSarbitrary-method
#'   wetted_area,CScircle-method
#' @description Calculates the wetted area of a CSarbitrary or CScircle
#'   object for given water levels.
#' @usage wetted_area(object, h, ret = "A")
#' @param object An object of class CSarbitrary or CScircle.
#' @param h A numeric vector of water levels (m). For CScircle, only a
#'   single numeric value is allowed.
#' @param ret A character string. If `A`, returns total wetted area.
#'   If `Aii`, returns wetted area by segment.
#' @return A numeric vector or matrix of wetted areas based on the `ret`
#'   argument.
#' @importFrom methods new validObject
#' @importFrom stats approx
#' @examples
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(x = x, z = z, xb_l = 4, xb_r = 9, kSt_B = 35,
#'                   kSt_l = 45, kSt_r = 45)
#'
#' # Calculate total wetted area at water levels 1 m and 2 m
#' h <- c(1, 2)
#' wetted_area(cs, h, ret = "A")
#'
#' # Calculate wetted area for each segment at the same water levels
#' wetted_area(cs, h, ret = "Aii")
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 1, kSt = 75)
#'
#' # Calculate total wetted area at water level 1 m
#' h <- 1
#' wetted_area(csC, h)
#' @export


setMethod("wetted_area", "CSarbitrary", function(object, h, ret = "A") {
  n_segments <- length(x(object)) - 1  # Number of segments

  # Initialize A_h depending on return type
  A_h <- if (ret == "A") rep(NA, length(h)) else matrix(nrow = length(h),
                                                        ncol = n_segments)

  for (jj in seq_along(h)) {
    zh <- min(z(object)) + h[jj]  # Level of h [m a.s.l.]
    Aii <- numeric(n_segments)    # Vector for areas


    for (ii in seq_len(n_segments)) {
      if (z(object)[ii] < zh & z(object)[ii + 1] < zh) {
        # Both points below water level
        Aii[ii] <- 0.5 * ((zh - z(object)[ii]) + (zh - z(object)[ii + 1])) *
          (x(object)[ii + 1] - x(object)[ii])
      } else if (z(object)[ii] >= zh & z(object)[ii + 1] < zh) {
        # One point below water level (riverbank)
        x_inters <- approx(x = z(object)[ii:(ii + 1)], y = x(object)[ii:(ii + 1)],
                           xout = zh)$y
        Aii[ii] <- 0.5 * (zh - z(object)[ii + 1]) * (x(object)[ii + 1] - x_inters)
      } else if (z(object)[ii] < zh & z(object)[ii + 1] >= zh) {
        # One point below water level (riverbank)
        x_inters <- approx(x = z(object)[ii:(ii + 1)], y = x(object)[ii:(ii + 1)],
                           xout = zh)$y
        Aii[ii] <- 0.5 * (zh - z(object)[ii]) * (x_inters - x(object)[ii])
      } else {
        # Both points above water level (dry)
        Aii[ii] <- 0
      }
    }

    if (ret == "A") {
      A_h[jj] <- sum(Aii)
    } else {
      A_h[jj, ] <- Aii
    }
  }

  return(A_h)
})


# Wetted Perimeter
#------------------------------------------------------------------------------
#' @title Wetted Perimeter
#' @name wetted_perimeter
#' @aliases wetted_perimeter wetted_perimeter,CSarbitrary-method
#'   wetted_perimeter,CScircle-method
#' @description Calculates the wetted perimeter of a CSarbitrary or CScircle
#'   object for given water levels.
#' @usage wetted_perimeter(object, h, ret = "P")
#' @param object An object of class CSarbitrary or CScircle.
#' @param h A numeric vector of water levels (m). For CScircle, only a
#'   single numeric value is allowed.
#' @param ret A character string. If `P`, returns total wetted perimeter.
#'   If `Pii`, returns wetted perimeter by segment.
#' @return A numeric vector or matrix of wetted perimeter based on the `ret`
#'   argument.
#' @importFrom methods new validObject
#' @importFrom stats approx
#' @examples
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(x = x, z = z, xb_l = 4, xb_r = 9, kSt_B = 35,
#'                   kSt_l = 45, kSt_r = 45)
#'
#' # Calculate total wetted perimeter at water levels 1 m and 2 m
#' h <- c(1, 2)
#' wetted_perimeter(cs, h, ret = "P")
#'
#' # Calculate wetted perimeter for each segment at the same water levels
#' wetted_perimeter(cs, h, ret = "Pii")
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 1, kSt = 75)
#'
#' # Calculate total wetted perimeter at water level 1 m
#' h <- 1
#' wetted_perimeter(csC, h)
#' @export


setMethod("wetted_perimeter", "CSarbitrary", function(object, h, ret = "P") {
  n_segments <- length(x(object)) - 1  # Number of segments

  # Initialize P_h depending on return type
  P_h <- if (ret == "P") rep(NA, length(h)) else matrix(nrow = length(h),
                                                        ncol = n_segments)

  for (jj in seq_along(h)) {
    zh <- min(z(object)) + h[jj]  # Level of h [m a.s.l.]
    Pii <- numeric(n_segments)    # Vector for hydraulic radius


    for (ii in seq_len(n_segments)) {
      if (z(object)[ii] < zh & z(object)[ii + 1] < zh) {
        # Both points below water level
        Pii[ii] <- sqrt((x(object)[ii + 1] - x(object)[ii])^2 +
                          (z(object)[ii + 1] - z(object)[ii])^2)
      } else if (z(object)[ii] >= zh & z(object)[ii + 1] < zh) {
        # One point below water level (riverbank)
        x_inters <- approx(x = z(object)[ii:(ii + 1)], y = x(object)[ii:(ii + 1)],
                           xout = zh)$y
        Pii[ii] <- sqrt((x(object)[ii + 1] - x_inters)^2 +
                          (z(object)[ii + 1] - zh)^2)
      } else if (z(object)[ii] < zh & z(object)[ii + 1] >= zh) {
        # One point below water level (riverbank)
        x_inters <- approx(x = z(object)[ii:(ii + 1)], y = x(object)[ii:(ii + 1)],
                           xout = zh)$y
        Pii[ii] <- sqrt((x_inters - x(object)[ii])^2 +
                          (zh - z(object)[ii])^2)
      } else {
        # Both points above water level (dry)
        Pii[ii] <- 0
      }
    }

    if (ret == "P") {
      P_h[jj] <- sum(Pii)
    } else {
      P_h[jj, ] <- Pii
    }
  }

  return(P_h)
})


# calculate mean roughness (Einstein)
#------------------------------------------------------------------------------
#' @title Mean Roughness
#' @name mean_roughness
#' @aliases mean_roughness,CSarbitrary-method
#' @description Calculates the mean roughness of a CSarbitrary object for a
#'  given
#'   set of water levels, based on Einstein (1934).
#' @usage mean_roughness(object, h)
#' @param object A CSarbitrary object.
#' @param h A numeric vector of water levels [m].
#' @return A numeric vector representing the mean roughness for the given water
#'  levels.
#' @examples
#' # Example usage:
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(x = x, z = z, xb_l = 4, xb_r = 9, kSt_B = 35,
#'                   kSt_l = 45, kSt_r = 45)
#' h_levels <- c(1, 2)  # water levels
#' mean_roughness(cs, h_levels)
#' @export
#'
setMethod(
  "mean_roughness", "CSarbitrary",
  function(object, h) {
    # Check for required numeric inputs
    if (!is.numeric(kSt_B(object))) {
      warning("kSt_B is missing")
      return(NULL)
    }
    if (!is.numeric(kSt_l(object))) {
      warning("kSt_l is missing")
      return(NULL)
    }
    if (!is.numeric(kSt_r(object))) {
      warning("kSt_r is missing")
      return(NULL)
    }

    # Calculate wetted perimeters
    Pii <- wetted_perimeter(object, h = h, ret = "sep")
    xm <- 0.5 * (
      x(object)[2:length(x(object))] + x(object)[1:(length(x(object)) - 1)]
    )

    # Separate wetted perimeters for different regions
    Pl <- rowSums(Pii[, xb_l(object) > xm, drop = FALSE])
    Ps <- rowSums(Pii[, xb_l(object) < xm & xb_r(object) > xm, drop = FALSE])
    Pr <- rowSums(Pii[, xb_r(object) < xm, drop = FALSE])

    # Total wetted perimeter
    P <- rowSums(Pii)

    # Mean roughness (Einstein)
    kSt_m <- P^(2 / 3) / (
      (Pl / kSt_l(object)^(3 / 2)) +
        (Ps / kSt_B(object)^(3 / 2)) +
        (Pr / kSt_r(object)^(3 / 2))
    )^(2 / 3)

    return(kSt_m)
  }
)

# calculate flow velocity
#------------------------------------------------------------------------------
#' @title Flow Velocity
#' @name flow_velocity
#' @aliases  flow_velocity flow_velocity,CSarbitrary-method flow_velocity,CScircle-method
#' @description Calculates the flow velocity of a CSarbitrary or CScircle object for a
#'  given water level and bottom slope under uniform flow conditions.
#' @usage flow_velocity(object, h, J, method = "Strickler",nu = 1.14e-6,...)
#' @param object A CSarbitrary or CScircle object.
#' @param h Flow depth [m].
#' @param J Bottom slope  [-].
#' @param method Method to calculate the roughness. Allowed are "Strickler"
#'  (equal roughness) "Einstein" (mean roughness) and "Prandtl-Coolebrook-White".
#' @param nu Kinematic viscosity [m2/s]. Only for CScircle objects
#' @param ... Additional arguments.
#' @return Flow velocity [m/s]
#' @examples
#'
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(x = x, z = z, xb_l = 4, xb_r = 9, kSt_B = 35,
#'                   kSt_l = 45, kSt_r = 45)
#' flow_velocity(cs, h = 1,J = 0.01, method = "Einstein")
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 0.7,ks = 1.5, kSt = 75)
#' flow_velocity(csC, h = 0.46, J = 0.004)
#' flow_velocity(csC, h = 0.46, J = 0.004, method = "Prandtl-Coolebrook-White")
#'
#' @export
#'


setMethod(
  "flow_velocity", "CSarbitrary",
  function(object, h, J, method = "Strickler",nu = 1.14e-6) {
    # Define variables for velocity calculation
    A <- wetted_area(object, h = h)
    P <- wetted_perimeter(object, h = h)

    # Equal roughness across the entire cross-section
    if (method == "Strickler") {
      if (!is.numeric(kSt_B(object))) {
        warning("kSt_B is missing")
        return(NULL)
      }
      v <- kSt_B(object) * sqrt(J) * (A / P)^(2 / 3)

      # Mean Strickler roughness (Einstein)
    } else if (method == "Einstein") {
      kSt_m <- mean_roughness(object, h = h)
      v <- kSt_m * sqrt(J) * (A / P)^(2 / 3)

    } else {
      stop("Unknown method. Use 'Strickler' or 'Einstein'.")
    }

    return(v)
  }
)


# Calculate coordinates of water level for plotting
#------------------------------------------------------------------------------
setMethod(
  "WL_coords", "CSarbitrary",
  function(object, h, v) {

    z_talweg <- min(z(object))
    z_wsp <- z_talweg + h
    pl <- list()
    pr <- list()

    # Left points
    zmax_l <- max(z(object)[1:which(x(object) %in% xb_l(object))])

    if (h > zmax_l) {
      pl$x <- zmax_l
      pl$y <- min(x(object))
    } else {
      xxx <- z_wsp > z(object)
      beg <- which(diff(xxx) == 1)[1]
      pl <- approx(
        x = z(object)[beg:(beg + 1)],
        y = x(object)[beg:(beg + 1)],
        xout = z_wsp
      )
    }

    # Right points
    zmax_r <- max(z(object)[which(x(object) %in% xb_r(object)):
                              length(z(object))])

    if (h > zmax_r) {
      pr$x <- zmax_r
      pr$y <- max(x(object))
    } else {
      xxx <- z_wsp > z(object)
      end <- tail(which(diff(xxx) == -1), 1)
      pr <- approx(
        x = z(object)[end:(end + 1)],
        y = x(object)[end:(end + 1)],
        xout = z_wsp
      )
    }

    pb <- c(
      pl$y, pl$y,
      x(object)[which(x(object) >= pl$y & x(object) <= pr$y)],
      pr$y, pr$y
    )

    pz <- c(
      z_wsp, pl$x,
      z(object)[which(x(object) >= pl$y & x(object) <= pr$y)],
      pr$x, z_wsp
    )

    return(list(zmax_l = zmax_l, zmax_r = zmax_r, pl = pl, pr = pr, pb = pb,
                pz = pz))
  }
)

# Calculate coordinates of energy level for plotting
#------------------------------------------------------------------------------
setMethod(
  "EL_coords", "CSarbitrary",
  function(object, h, v) {

    z_talweg <- min(z(object))
    z_wsp <- z_talweg + h
    z_EL <- z_wsp + v^2 / (2 * 9.81)
    pl <- list()
    pr <- list()

    # Left points
    zmax_l <- max(z(object)[1:which(x(object) %in% xb_l(object))])

    if (z_EL > zmax_l) {
      pl$x <- zmax_l
      pl$y <- min(x(object))
    } else {
      xxx <- z_EL > z(object)
      beg <- which(diff(xxx) == 1)[1]
      pl <- approx(
        x = z(object)[beg:(beg + 1)],
        y = x(object)[beg:(beg + 1)],
        xout = z_EL
      )
    }

    # Right points
    zmax_r <- max(z(object)[which(x(object) %in%
                                    xb_r(object)):length(z(object))])

    if (z_EL > zmax_r) {
      pr$x <- zmax_r
      pr$y <- max(x(object))
    } else {
      xxx <- z_EL > z(object)
      end <- tail(which(diff(xxx) == -1), 1)
      pr <- approx(
        x = z(object)[end:(end + 1)],
        y = x(object)[end:(end + 1)],
        xout = z_EL
      )
    }

    return(list(pl = pl, pr = pr))
  }
)


# Froude number (Fr)
#------------------------------------------------------------------------------
#' @title Froude Number
#' @name froude_number
#' @aliases  froude_number froude_number,CSarbitrary-method froude_number,CScircle-method
#' @description Calculates the froude number of a CSarbitrary or CScircle object for a
#'  given water level and velocity under uniform flow conditions.
#' @usage froude_number(object, v, h)
#' @param object A CSarbitrary or CScircle object.
#' @param h Flow depth [m].
#' @param v Flow velocity  [m/s].
#' @return Froude number [-]
#' @examples
#'
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(x = x, z = z, xb_l = 4, xb_r = 9, kSt_B = 35,
#'                   kSt_l = 45, kSt_r = 45)
#' froude_number(cs,h=1, v = 2.5)
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 0.7,ks = 1.5, kSt = 75)
#' froude_number(csC, h = 0.46, v = 2.5)
#'
#' @export
#'
setMethod(
  "froude_number", "CSarbitrary",
  function(object, v, h) {

    A <- wetted_area(object, h = h)
    delta_A <- (wetted_area(object, h = h + 0.01) -
                  wetted_area(object, h = h - 0.01)) / (2 * 0.01)

    Fr <- v/sqrt(9.81*A/delta_A)

    return(Fr)
  }
)




# Calculate Flow Depth
#------------------------------------------------------------------------------
#' @title Flow Depth
#' @name flow_depth
#' @aliases flow_depth flow_depth,CSarbitrary-method flow_depth,CScircle-method
#' @description Calculates the flow depth of a CSarbitrary or CScircle object for
#'   a given discharge and bottom slope under uniform flow conditions.
#' @usage flow_depth(object, Q, J, method = "Strickler", ret = "all", plot = FALSE)
#' @param object A CSarbitrary or CScircle object.
#' @param Q Discharge [m3/s].
#' @param J Bottom slope [-].
#' @param method Method to calculate the roughness. Allowed are "Strickler"
#'  (equal roughness) "Einstein" (mean roughness) and "Prandtl-Coolebrook-White".
#' @param ret Defines the result returned by the function.
#' @param plot Logical; if `TRUE`, plots the results.
#' @return A list containing the following hydraulic variables:
#' \describe{
#'   \item{h}{Flow depth [m].}
#'   \item{v}{Flow velocity [m/s].}
#'   \item{Fr}{Froude number [-].}
#'   \item{kSt_m}{Mean roughness [m^(1/3)/s] (if method = "Einstein").}
#'   \item{A}{Wetted area [m^2].}
#'   \item{P}{Wetted perimeter [m].}
#' }
#' @examples
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(
#'   x = x, z = z, xb_l = 4, xb_r = 9,
#'   kSt_B = 35, kSt_l = 45, kSt_r = 45
#' )
#' flow_depth(cs, Q = 8.677, J = 0.0001, method = "Einstein", ret = "h")
#' flow_depth(cs, Q = 8.677, J = 0.0001, method = "Einstein", plot = TRUE)
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 0.7, ks = 1.5, kSt = 75)
#' flow_depth(csC, Q = 0.46, J = 0.004)
#' flow_depth(csC, Q = 0.46, J = 0.004, method = "Prandtl-Coolebrook-White", plot = TRUE)
#'
#' @importFrom grDevices dev.off rgb
#' @importFrom graphics legend lines points polygon axis symbols
#' @importFrom stats uniroot
#' @importFrom utils tail
#' @export


setMethod(
  "flow_depth", "CSarbitrary",
  function(object, Q, J, method = "Strickler", ret = "all", plot = FALSE) {

    # Define function to find root
    fcn <- function(h) {
      A <- wetted_area(object, h = h)
      P <- wetted_perimeter(object, h = h)
      v <- flow_velocity(object, h = h, J = J, method = method)
      Q / A - v
    }

    # Find root of the function
    h <- uniroot(fcn, interval = c(1e-5, 10000), check.conv = TRUE)$root
    v <- flow_velocity(object, h = h, J = J, method = method)

    if(h>max(z(object))){
      warning("h exceeds dike level")
    }

    # Plot if required
    if (plot) {
      try(dev.off(), silent = TRUE)

      z_min <- min(z(object))
      z_wsp <- z_min + h
      z_el <- z_wsp + v^2 / (2 * 9.81)
      ylim<-if((min(z(object))+h+ v^2/(2*9.81))<max(z(object))){
        c(0,max(z(object)))
      } else{c(0,min(z(object))+h+ v^2/(2*9.81))}

      # Plot cross-section and water profile
      plot(
        x(object), z(object), type = "l", xlab = "x [m]", ylab = "z [m]",
        main = "", cex.main = 0.8, asp = 1, ylim = ylim,)
      lines(
        x = c(
          WL_coords(object, h = h, v = v)$pl$y,
          WL_coords(object, h = h, v = v)$pr$y
        ),
        y = rep(z_wsp, 2), col = "blue"
      )
      polygon(
        WL_coords(object, h = h, v = v)$pb,
        WL_coords(object, h = h, v = v)$pz,
        col = rgb(0.68, 0.85, 0.9, 0.3), border = NA
      )
      lines(
        x = c(
          EL_coords(object, h = h, v = v)$pl$y,
          EL_coords(object, h = h, v = v)$pr$y
        ),
        y = rep(z_el, 2), col = "red", lty = 4
      )
      legend(
        "bottomleft", inset = c(0, 1), xpd = TRUE,
        legend = c(
          "CS", "water level", "energy line", "bank bottom",
          paste("Q =", Q, "m3/s"), paste("J =", J),
          paste("v =", round(v, 2), "m/s"),
          paste("F =", round(froude_number(object, v = v, h = h), 2),
                ", h =", round(h, 2), "m")
        ),
        lty = c(1, 1, 4, NA, NA, NA, NA, NA),
        pch = c(NA, NA, NA, 6, NA, NA),
        col = c("black", "blue", "red", "black", NA, NA, NA),
        bty = "n", cex = 0.8, ncol = 2
      )
      points(
        c(xb_l(object), xb_r(object)),
        c(
          z(object)[x(object) %in% xb_l(object)],
          z(object)[x(object) %in% xb_r(object)]
        ),
        pch = 6
      )

    }

    # Output results based on the ret argument
    if (ret == "all") {
      res <- list(
        h = h,
        v = v,
        Fr = froude_number(object, v = v, h = h),
        A = wetted_area(object, h = h),
        P = wetted_perimeter(object, h = h)
      )
      if (method == "Einstein") {
        res$kSt_m <- mean_roughness(object, h = h)
      }
      return(res)
    } else if (ret == "h") {
      return(h)
    } else if (ret == "v") {
      return(v)
    }
  }
)

# Flow
#------------------------------------------------------------------------------
#' @title Flow
#' @name flow
#' @aliases flow flow,CSarbitrary-method flow,CScircle-method
#' @description Calculates the discharge of a CSarbitrary or CScircle object for
#'   a given flow depth and bottom slope under uniform flow conditions.
#' @usage flow(object, h, J, method = "Strickler", ret = "all", plot = FALSE)
#' @param object A CSarbitrary or CScircle object.
#' @param h Flow depth [m].
#' @param J Bottom slope [-].
#' @param method Method to calculate the roughness. Allowed are "Strickler"
#'  (equal roughness) "Einstein" (mean roughness) and "Prandtl-Coolebrook-White".
#' @param ret Defines the result returned by the function.
#' @param plot Logical; if `TRUE`, plots the results.
#' @return A list containing the following hydraulic variables:
#' \describe{
#'   \item{Q}{Discharge [m3/s].}
#'   \item{v}{Flow velocity [m/s].}
#'   \item{kSt_m}{Mean roughness [m^(1/3)/s] (if method = "Einstein").}
#'   \item{A}{Wetted area [m^2].}
#' }
#' @examples
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(
#'   x = x, z = z, xb_l = 4, xb_r = 9,
#'   kSt_B = 35, kSt_l = 45, kSt_r = 45
#' )
#' flow(cs, h = 2, J = 0.0001, method = "Einstein", ret = "Q")
#' flow(cs, h = 2, J = 0.0001, method = "Einstein", plot = TRUE)
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 0.7, ks = 1.5, kSt = 75)
#' flow(csC, h = 0.46, J = 0.004)
#' flow(csC, h = 0.46, J = 0.004, method = "Prandtl-Coolebrook-White", plot = TRUE)
#' @export
#'
setMethod(
  "flow", "CSarbitrary",
  function(object, h, J, method = "Strickler", ret = "all", plot = FALSE) {

    if(h>max(z(object))){
      warning("h exceeds dike level")
    }

    A <- wetted_area(object, h = h)
    v <- flow_velocity(object, h = h, J = J, method = method)
    Q <- v * A

    # Plot section
    if (plot) {
      try(dev.off(), silent = TRUE)

      z_min <- min(z(object))
      ylim <- if ((z_min + h + v^2 / (2 * 9.81)) < max(z(object))) {
        c(0, max(z(object)))
      } else {
        c(0, z_min + h + v^2 / (2 * 9.81))
      }

      plot(
        x(object), z(object), type = "l", xlab = "x [m]", ylab = "z [m]",
        main = "", cex.main = 0.8, asp = 1, ylim = ylim
      )

      # Water table
      z_wsp <- z_min + h
      wl <- WL_coords(object, h = h, v = v)
      lines(
        x = c(wl$pl$y, wl$pr$y), y = rep(z_wsp, 2), col = "blue"
      )
      polygon(
        wl$pb, wl$pz, col = rgb(0.68, 0.85, 0.9, 0.3), border = NA
      )

      # Energy line
      z_el <- z_wsp + v^2 / (2 * 9.81)
      el <- EL_coords(object, h = h, v = v)
      lines(
        x = c(el$pl$y, el$pr$y), y = rep(z_el, 2), col = "red", lty = 4
      )

      legend(
        "bottomleft", inset = c(0, 1), xpd = TRUE,
        legend = c(
          "CS", "water level", "energy line", "bank bottom",
          paste("Q =", round(Q, 2), "m^3/s"),
          paste("J =", J),
          paste("v =", round(v, 2), "m/s"),
          paste("F =", round(froude_number(object, v = v, h = h), 2))
        ),
        lty = c(1, 1, 2, NA, NA, NA, NA, NA), pch = c(NA, NA, NA, 6, NA, NA),
        col = c("black", "blue", "red", "black", NA, NA), bty = "n",
        cex = 0.8, ncol = 2
      )

      points(
        c(xb_l(object), xb_r(object)),
        c(z(object)[x(object) %in% xb_l(object)],
          z(object)[x(object) %in% xb_r(object)]),
        pch = 6
      )

    }

    # Output section
    if (ret == "all") {
      res <- list(Q = Q, v = v, A = A)
      if (method == "Einstein") {
        res$kSt_m <- mean_roughness(object, h = h)
      }
      return(res)
    } else if (ret == "Q") {
      return(Q)
    } else if (ret == "v") {
      return(v)
    }
  }
)



# Flow_max
#------------------------------------------------------------------------------
#' @title Maximum Flow
#' @name flow_max
#' @aliases flow_max flow_max,CSarbitrary-method flow_max,CScircle-method
#' @description Calculates the maximum discharge of a CSarbitrary or CScircle
#' object for a given bottom slope under uniform flow conditions.
#' @usage flow_max(object, J, method = "Strickler", ret = "all", plot = FALSE)
#' @param object A CSarbitrary or CScircle object.
#' @param J Bottom slope [-].
#' @param method Method to calculate the roughness. Allowed are "Strickler"
#'  (equal roughness) "Einstein" (mean roughness) and "Prandtl-Coolebrook-White".
#' @param ret Defines the result returned by the function.
#' @param plot Logical; if TRUE, plots the results.
#' @return A list containing the following hydraulic variables:
#' \describe{
#'   \item{Qmax}{Maximum discharge [m3/s].}
#'   \item{hmax}{Maximum flow depth [m].}
#'   \item{v}{Flow velocity [m/s].}
#'   \item{kSt_m}{Mean roughness [m^(1/3)/s] (if method = "Einstein").}
#'   \item{A}{Wetted area [m2].}
#' }
#' @examples
#' # Example for CSarbitrary object
#' x <- c(0, 4, 9, 13)
#' z <- c(2, 0, 0, 2)
#' cs <- CSarbitrary(
#'   x = x, z = z, xb_l = 4, xb_r = 9,
#'   kSt_B = 35, kSt_l = 45, kSt_r = 45
#' )
#' flow_max(cs, J=0.0001, method="Einstein",ret="Qmax")
#' flow_max(cs, J=0.0001, method="Einstein",plot=TRUE)
#'
#' # Example for CScircle object
#' csC <- CScircle(Di = 0.7, ks = 1.5, kSt = 75)
#' flow_max(csC, J=0.004)
#' flow_max(csC, J = 0.004, method = "Prandtl-Coolebrook-White", plot = TRUE)
#' @export
#'


setMethod(
  "flow_max", "CSarbitrary",
  function(object, J, method = "Strickler", ret = "all", plot = FALSE) {

    # Check if xb_l or xb_r is missing
    if (!is.numeric(xb_l(object))) {
      warning("xb_l is missing")
      return(NULL)
    }
    if (!is.numeric(xb_r(object))) {
      warning("xb_r is missing")
      return(NULL)
    }

    # Calculate maximal possible level
    zmax_l <- max(z(object)[xb_l(object) >= x(object)], na.rm = TRUE)
    zmax_r <- max(z(object)[xb_r(object) <= x(object)], na.rm = TRUE)
    zmax <- min(zmax_l, zmax_r)
    hmax <- zmax - min(z(object), na.rm = TRUE)

    # Calculate flow and velocity for maximal level
    Qmax <- flow(object, hmax, J = J, method = method, ret = "Q")
    v <- flow_velocity(object, h = hmax, J = J, method = method)

    # Plotting
    if (plot==TRUE) {
      suppressWarnings(try(dev.off(), silent = TRUE))

      ylim <- if ((min(z(object)) + hmax + v^2 / (2 * 9.81)) < max(z(object))) {
        c(0, max(z(object)))
      } else {
        c(0, min(z(object)) + hmax + v^2 / (2 * 9.81))
      }

      plot(
        x = x(object), y = z(object), type = "l",
        xlab = "x [m]", ylab = "z [m]",
        main = "Qmax", cex.main = 0.8, asp = 1,
        ylim = ylim
        )

      # Water level and energy line
      z_talweg <- min(z(object), na.rm = TRUE)
      z_wsp <- z_talweg + hmax
      lines(
        x = c(
          WL_coords(object, h = hmax, v = v)$pl$y,
          WL_coords(object, h = hmax, v = v)$pr$y
        ),
        y = rep(z_wsp, 2), col = "blue"
      )
      polygon(
        x = WL_coords(object, h = hmax, v = v)$pb,
        y = WL_coords(object, h = hmax, v = v)$pz,
        col = rgb(0.68, 0.85, 0.9, 0.3), border = NA
      )
      z_el <- z_wsp + v^2 / (2 * 9.81)
      lines(
        x = c(
          EL_coords(object, h = hmax, v = v)$pl$y,
          EL_coords(object, h = hmax, v = v)$pr$y
        ),
        y = rep(z_el, 2), col = "red", lty = 4
      )

      legend(
        "bottomleft", inset = c(0, 1), xpd = TRUE,
        legend = c(
          "CS", "water level", "energy line", "bank bottom",
          paste("Q =", round(Qmax, 2), "m3/s"),
          paste("J =", J),
          paste("v =", round(v, 2), "m/s"),
          paste("F =", round(froude_number(object, v = v, h = hmax), 2))
        ),
        lty = c(1, 1, 2, NA, NA, NA, NA, NA), pch = c(NA, NA, NA, 6, NA, NA),
        col = c("black", "blue", "red", "black", NA, NA, NA),
        bty = "n", cex = 0.8, ncol = 2
      )

      points(
        x = c(xb_l(object), xb_r(object)),
        y = c(z(object)[x(object) %in% xb_l(object)],
              z(object)[x(object) %in% xb_r(object)]),
        pch = 6
      )
    }
    # Return values
    switch(
      ret,
      "all" = list(
        Qmax = Qmax, hmax = hmax, v = v,
        A = wetted_area(object, h = hmax),
        kSt_m = if (method == "Einstein") mean_roughness(object, h = hmax)
      ),
      "Qmax" = Qmax,
      "hmax" = hmax,
      "v" = flow_velocity(object, h = hmax, J = J, method = method)
    )
  }
)

#------------------------------------------------------------------------------
# flow_max_freeboard
#------------------------------------------------------------------------------
#' @title Maximum Flow Including Freeboard
#' @name flow_max_freeboard
#' @description Calculates the maximum discharge of a CSarbitrary
#' object including a freebord for a given bottom slope under uniform flow conditions.
#' @aliases flow_max_freeboard,CSarbitrary-method
#' @usage flow_max_freeboard(object, J, type = "KOHS", sigma_wz = 0, fw = TRUE, fv = FALSE, ft = 0,
#' fe = NULL, fe_min = 0, fe_max = Inf, method = "Strickler",
#' ret = "all", plot = FALSE)
#' @param object A CSarbitrary object.
#' @param J Bottom slope [-].
#' @param type Type of freeboard calculation. Defaults to "KOHS".
#' @param sigma_wz Uncertainty in bed elevation (morphodynamics) [m].
#' @param fw Logical; considers freeboard due to uncertainty in water elevation.
#'   If TRUE, calculates according to KOHS; if FALSE, sets fw = 0.
#' @param fv Logical; considers freeboard due to waves. If `TRUE`, calculates
#'   according to KOHS; if FALSE, sets fv = 0.
#' @param ft Freeboard due to driftwood based on KOHS (2013) [m].
#' @param fe Fixed freeboard value to override calculations [m].
#' @param fe_min Minimum freeboard [m].
#' @param fe_max Maximum freeboard [m].
#' @param method Method to calculate the roughness. Allowed are "Strickler"
#'  (equal roughness) and "Einstein" (mean roughness).
#' @param ret Definition of the result returned by the function ("all", "Qmax",
#'   "hmax", "fe", or "v").
#' @param plot Logical; whether to plot the results.
#' @return Depending on ret, returns flow, water level, velocity, or all
#'   details.
#' @references KOHS (2013). Freibord bei Hochwasserschutzprojekten und
#'   Gefahrenbeurteilungen - Empfehlungen der Kommission Hochwasserschutz KOHS.
#'   Wasser Energie Luft 105(1): 43-53.
#' @examples
#' # Cross section
#' x <- c(-0.85, 3, 15, 18.85)
#' z <- c(3.85, 0, 0, 3.85)
#' cs<- CSarbitrary(x = x, z = z, xb_l = 3, xb_r = 15,
#'                                       kSt_B = 45)
#'
#' # Channel
#' flow_max_freeboard(cs, sigma_wz = 0.3, fv = FALSE, J = 2.2 * 10^-2)
#' # Dam
#' flow_max_freeboard(cs, sigma_wz = 0.3, fv = TRUE, J = 2.2 * 10^-2)
#' # Bridge
#' flow_max_freeboard(cs, sigma_wz = 0.3, fv = TRUE, ft = 0.5,
#'            J = 2.2 * 10^-2)
#'
#' # Sensitivity analysis for slope
#' J <- seq(1, 3, 0.1) * 10^-2
#' Q <- sapply(J, function(J) {
#'   flow_max_freeboard(cs, sigma_wz = 0.3, fv = TRUE, ft = 0.5,
#'              J = J)$Qmax
#' })
#' plot(J, Q, type = "l")
#' @export
#'
#'
setMethod(
  "flow_max_freeboard",
  "CSarbitrary",
  function(object, J, type = "KOHS", sigma_wz = 0, fw = TRUE, fv = FALSE, ft = 0,
           fe = NULL, fe_min = 0, fe_max = Inf, method = "Strickler",
           ret = "all", plot = FALSE) {

    # Check and validate inputs
    if (!is.numeric(xb_l(object))) {
      warning("Left bank xb_l is missing.")
      return(NULL)
    }
    if (!is.numeric(xb_r(object))) {
      warning("Right bank xb_r is missing.")
      return(NULL)
    }

    # Get maximal levee levels
    zmax_l <- max(z(object)[xb_l(object) >= x(object)], na.rm = TRUE)
    zmax_r <- max(z(object)[xb_r(object) <= x(object)], na.rm = TRUE)
    zmax <- min(zmax_l, zmax_r, na.rm = TRUE)

    # Constant freeboard
    if (type == "const_fe") {
      if (!is.numeric(fe)) {
        warning("Fixed freeboard 'fe' is not defined.")
        return(NULL)
      }
      hmax <- zmax - fe - min(z(object))

    } else if (type == "KOHS") {
      # Without bridge
      calc_freeboard <- function(h, zmax) {
        v <- flow_velocity(object, h = h, J = J, method = method)
        fe <- freeboard(v, h, sigma_wz, fw, fv, ft, min = fe_min, max = fe_max)
        hmax <- zmax - fe - min(z(object))
        h - hmax
      }
      hmax_levee <- uniroot(
        calc_freeboard, interval = c(10^-5, 1e4), zmax = zmax, check.conv = TRUE
      )$root
      fe_levee <- zmax - min(z(object)) - hmax_levee

      hmax <- hmax_levee
      fe <- fe_levee
      fe_type <- "levee"


    } else {
      warning("Unknown 'type' parameter.")
      return(NULL)
    }

    # Calculate flow and velocity
    Qmax <- flow(object, hmax, J = J, method = method, ret = "Q")
    v <- flow_velocity(object, h = hmax, J = J, method = method)

    # Optional plot
    if (plot) {
      plot(x(object), z(object), type = "l", xlab = "x [m]", ylab = "z [m]",
           main = "Qfreeboard", ylim = c(0, max(c(z(object), hmax + v^2 / 19.62))))
      lines(x = c(WL_coords(object, h = hmax, v = v)$pl$y,
                  WL_coords(object, h = hmax, v = v)$pr$y),
            y = rep(min(z(object)) + hmax, 2), col = "blue")
        }

    # Prepare outputs
    if (ret == "all") {
      result <- list(
        Qmax = Qmax,
        hmax = hmax,
        fe = fe,
        v = v,
        A = wetted_area(object, h = hmax)
      )
      if (method == "Einstein") {
        result$kSt_m <- mean_roughness(object, h = hmax)
      }
      if (type == "KOHS") {
        result$fe_type <- fe_type
      }
      return(result)
    }
    switch(ret, Qmax = Qmax, hmax = hmax, fe = fe, v = v)
  }
)

