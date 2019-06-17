
#' Add titer impulses to a 3d landscape plot
#'
#' @param data3js The r3js data object
#' @param object The antibody landscape object
#' @param zlim z limits for the plot
#' @param ... Additional plotting parameters
#' 
#' @export
#' 
lndscp3d_titres <- function(data3js,
                            object,
                            zlim,
                            ...) {
  
  # Get the plot pars
  pars <- ablandscape.par(...)
  
  # Define NA titres
  na_titers <- is.na(object$logtiters)
  
  # Get the coordinates
  x <- object$coords[,1]
  y <- object$coords[,2]
  z <- object$logtiters
  
  # Get the indices of the points fitted
  indices  <- object$ag_indices
  ag_cols  <- object$acmap$ag_cols_fill[indices]
  ag_names <- object$acmap$ag_names[indices]
  
  # Plot the impulses and points if requested.
  for (n in which(!na_titers)) {
    
    # Start a new group
    groupIDs <- c()
    
    # Plot basepoints
    data3js <- r3js::points3js(data3js,
                               x          = x[n],
                               y          = y[n],
                               z          = zlim[1]+0.01,
                               size       = pars$cex.titer,
                               col        = ag_cols[n],
                               highlight  = list(col = "red"),
                               label      = ag_names[n],
                               toggle     = "Titers",
                               depthWrite = FALSE,
                               dimensions = 2)
    groupIDs <- c(groupIDs, r3js::lastID(data3js))
    
    data3js <- r3js::points3js(data3js,
                               x                   = x[n],
                               y                   = y[n],
                               z                   = zlim[1]+0.01,
                               size                = pars$cex.titer,
                               pch                 = 1,
                               col                 = "black",
                               highlight           = list(col = "red"),
                               toggle              = "Titers",
                               label               = ag_names[n],
                               polygonOffset       = TRUE,
                               polygonOffsetFactor = -1.0,
                               polygonOffsetUnits  = -1.0,
                               dimensions          = 2)
    groupIDs <- c(groupIDs, r3js::lastID(data3js))
    
    
    # Plot impulses
    data3js <- r3js::lines3js(data3js,
                              x = rep(x[n], 2),
                              y = rep(y[n], 2),
                              z = c(zlim[1], z[n]),
                              col = "grey20",
                              highlight = list(col = "red"),
                              interactive = FALSE,
                              toggle = "Titers",
                              # geometry = TRUE,
                              lwd = 4)
    groupIDs <- c(groupIDs, r3js::lastID(data3js))
    
    data3js <- r3js::points3js(data3js,
                               x         = x[n],
                               y         = y[n],
                               z         = z[n],
                               size      = pars$cex.titer,
                               col       = "grey20",
                               highlight = list(col = "red"),
                               label     = ag_names[n],
                               opacity   = 1,
                               toggle    = "Titers")
    groupIDs <- c(groupIDs, r3js::lastID(data3js))
    
    # Link the group
    data3js <- r3js::group3js(data3js, groupIDs)
    
  }
  
  invisible(data3js)
  
}
