
#' Generate a 3D map plot for an antibody landscape
#'
#' @param object The antibody landscape object
#' @param ... Graphical parameters
#'
#' @export
#'
lndscp3d_map <- function(object,
                         xlim,
                         ylim,
                         zlim,
                         ...) {
  
  # Get the plot pars
  pars <- ablandscape.par(...)
  
  # Open a new r3js plot
  data3js <- r3js::plot3js.new()
  
  ## Set plot window
  data3js <- r3js::plot3js.window(data3js,
                                  xlim = xlim,
                                  ylim = ylim,
                                  zlim = zlim,
                                  aspect = c(1, 1, 1))
  
  ## Set box
  data3js <- r3js::box3js(data3js,
                          type  = "l",
                          col   = "grey80",
                          sides = "xy")
  
  ## Display axes
  axis_labels <- zlim[1]:zlim[2]
  if(pars$zaxt == "linear"){
    axis_labels <- 2^axis_labels*10
    axis_labels[axis_labels == 5] <- "<10"
  }
  if(pars$zaxt == "log"){
    axis_labels[axis_labels == -1] <- "nd"
  }
  
  ## Add a base grid
  for(x in seq(from = min(xlim), to = max(xlim))){
    data3js <- r3js::lines3js(data3js,
                              x = c(x, x),
                              y = range(ylim),
                              z = c(zlim[1], zlim[1]),
                              col = pars$col.grid,
                              xpd = TRUE)
  }
  for(y in seq(from = floor(min(ylim)), to = ceiling(max(ylim)))){
    data3js <- r3js::lines3js(data3js,
                              x = range(xlim),
                              y = c(y, y),
                              z = c(zlim[1], zlim[1]),
                              col = pars$col.grid,
                              xpd = TRUE)
  }
  
  ## Add antigens and sera
  for (n in seq_len(nrow(object$acmap$ag_coords))) {
    data3js <- r3js::points3js(data3js,
                               x          = object$acmap$ag_coords[n,1],
                               y          = object$acmap$ag_coords[n,2],
                               z          = zlim[1]+0.005,
                               size       = pars$cex.basemap,
                               col        = object$acmap$ag_cols_fill[n],
                               highlight  = list(col = "red"),
                               label      = object$acmap$ag_names[n],
                               toggle     = "Basepoints",
                               depthWrite = FALSE,
                               dimensions = 2)
  }
  
  for (n in seq_len(nrow(object$acmap$sr_coords))) {
    data3js <- r3js::points3js(data3js,
                               x          = object$acmap$sr_coords[n,1],
                               y          = object$acmap$sr_coords[n,2],
                               z          = zlim[1]+0.005,
                               pch        = 0,
                               size       = pars$cex.basemap,
                               col        = object$acmap$sr_cols_outline[n],
                               highlight  = list(col = "red"),
                               label      = object$acmap$sr_names[n],
                               toggle     = "Basepoints",
                               depthWrite = FALSE,
                               dimensions = 2)
  }
  
  # Return new plot data
  invisible(data3js)
  
}
