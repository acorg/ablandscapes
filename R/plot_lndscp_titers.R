
#' Add titer impulses to a 3d landscape plot
#'
#' @param data3js The r3js data object
#' @param object The antibody landscape object
#' @param zlim z limits for the plot
#' @param ... Additional plotting parameters
#' 
#' @export
#' 
lndscp3d_titers <- function(
  data3js,
  object,
  zlim,
  show.impulses = TRUE,
  toggle = "Titers",
  titer_proportion_cutoff = 0.9,
  options = list()
  ) {
  
  # Get the plot pars
  pars <- do.call(ablandscape.par, options)
  
  # Get the coordinates
  x <- object$coords[,1]
  y <- object$coords[,2]
  if (is.null(dim(object$logtiters))) {
    z <- object$logtiters
  } else {
    z <- apply(object$titers, 2, \(x) meantiter::mean_titers(x, "truncated_normal", dilution_stepsize = 1)$mean)
    z_ntitrations <- colSums(!is.na(object$logtiters))
    z_excluded <- z_ntitrations / max(z_ntitrations) < titer_proportion_cutoff
    z[z_excluded] <- NA
    z[z < zlim[1]] <- zlim[1]
  }
  
  # Define NA titers
  na_titers <- is.na(z)
  
  # Take an average if z is a matrix
  if (!is.null(dim(z))) z <- colMeans(z)
  
  # Get the indices of the points fitted
  if (!is.null(object$acmap)) {
    indices  <- object$acmap_indices
    ag_cols  <- Racmacs::agFill(object$acmap)[indices]
    ag_names <- Racmacs::agNames(object$acmap)[indices]
  } else {
    ag_cols <- rep("grey50", length(z))
    if (!is.null(rownames(object$coords))) ag_names <- rownames(object$coords)
    else                                   ag_names <- paste("AG", seq_along(z))
  }
  
  # Plot the impulses and points if requested.
  for (n in which(!na_titers)) {
    
    # Start a new group
    groupIDs <- c()
    
    # Plot basepoints
    data3js <- r3js::points3js(data3js,
                               x          = x[n],
                               y          = y[n],
                               z          = zlim[1]+0.02,
                               size       = pars$cex.titer*2,
                               col        = ag_cols[n],
                               highlight  = list(col = "red"),
                               label      = ag_names[n],
                               toggle     = toggle,
                               depthWrite = FALSE,
                               dimensions = 2)
    groupIDs <- c(groupIDs, r3js::lastID(data3js))
    
    data3js <- r3js::points3js(data3js,
                               x                   = x[n],
                               y                   = y[n],
                               z                   = zlim[1]+0.02,
                               size                = pars$cex.titer*2,
                               pch                 = 1,
                               col                 = "black",
                               highlight           = list(col = "red"),
                               toggle              = toggle,
                               label               = ag_names[n],
                               polygonOffset       = TRUE,
                               polygonOffsetFactor = -1.0,
                               polygonOffsetUnits  = -1.0,
                               dimensions          = 2,
                               lwd                 = pars$lwd.titer.outline)
    groupIDs <- c(groupIDs, r3js::lastID(data3js))
    
    
    # Plot impulses
    if(show.impulses){
      data3js <- r3js::lines3js(data3js,
                                x = rep(x[n], 2),
                                y = rep(y[n], 2),
                                z = c(zlim[1], z[n]),
                                col = pars$col.impulse,
                                highlight = list(col = "red"),
                                interactive = FALSE,
                                toggle = toggle,
                                geometry = TRUE,
                                lwd = pars$lwd.impulse)
      groupIDs <- c(groupIDs, r3js::lastID(data3js))
      
      data3js <- r3js::points3js(data3js,
                                 x         = x[n],
                                 y         = y[n],
                                 z         = z[n],
                                 size      = pars$cex.titer.impulse,
                                 col       = "grey20",
                                 highlight = list(col = "red"),
                                 label     = ag_names[n],
                                 opacity   = 1,
                                 toggle    = toggle)
      groupIDs <- c(groupIDs, r3js::lastID(data3js))
    }
    
    # Link the group
    data3js <- r3js::group3js(data3js, groupIDs)
    
  }
  
  invisible(data3js)
  
}
