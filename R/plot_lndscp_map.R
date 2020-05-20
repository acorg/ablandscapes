
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
                         show.map.antigens = TRUE,
                         show.map.sera = TRUE,
                         aspect.z = 1,
                         toggle = "Basepoints",
                         ...) {
  
  # Get the map data
  map_ag_coords  <- Racmacs::agCoords(object$acmap)
  map_sr_coords  <- Racmacs::srCoords(object$acmap)
  map_ag_fill    <- Racmacs::agFill(object$acmap)
  map_sr_outline <- Racmacs::srOutline(object$acmap)
  map_ag_names   <- Racmacs::agNames(object$acmap)
  map_sr_names   <- Racmacs::srNames(object$acmap)
  
  # Get the plot pars
  pars <- ablandscape.par(...)
  
  # Open a new r3js plot
  data3js <- r3js::plot3js.new()
  
  ## Set plot window
  data3js <- r3js::plot3js.window(data3js,
                                  xlim = xlim,
                                  ylim = ylim,
                                  zlim = zlim,
                                  aspect = c(1, 1, aspect.z))
  
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
  if(show.map.antigens){
    for (n in seq_len(nrow(map_ag_coords))) {
      data3js <- r3js::points3js(data3js,
                                 x          = map_ag_coords[n,1],
                                 y          = map_ag_coords[n,2],
                                 z          = zlim[1]+0.005,
                                 size       = pars$cex.basemap*6,
                                 col        = map_ag_fill[n],
                                 opacity    = pars$opacity.basemap,
                                 highlight  = list(col = "red"),
                                 label      = map_ag_names[n],
                                 toggle     = toggle,
                                 depthWrite = FALSE,
                                 dimensions = 2)
    }
  }
  
  if(show.map.sera){
    for (n in seq_len(nrow(map_sr_coords))) {
      data3js <- r3js::points3js(data3js,
                                 x          = map_sr_coords[n,1],
                                 y          = map_sr_coords[n,2],
                                 z          = zlim[1]+0.005,
                                 pch        = 0,
                                 size       = pars$cex.basemap*6,
                                 col        = map_sr_outline[n],
                                 opacity    = pars$opacity.basemap,
                                 highlight  = list(col = "red"),
                                 label      = map_sr_names[n],
                                 toggle     = toggle,
                                 depthWrite = FALSE,
                                 dimensions = 2)
    }
  }
  
  # Return new plot data
  invisible(data3js)
  
}
