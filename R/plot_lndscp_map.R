
#' Generate a 3D map plot for an antibody landscape
#'
#' @param object The antibody landscape object
#' @param ... Graphical parameters
#'
#' @export
#'
lndscp3d_map <- function(
  data3js,
  fit,
  xlim,
  ylim,
  zlim,
  show.map.antigens = TRUE,
  show.map.sera = TRUE,
  aspect.z = 1,
  toggle = "Basepoints",
  options = list()
  ) {
  
  # Get the map data
  map_ag_coords  <- Racmacs::agCoords(fit$acmap)
  map_sr_coords  <- Racmacs::srCoords(fit$acmap)
  map_ag_fill    <- Racmacs::agFill(fit$acmap)
  map_sr_outline <- Racmacs::srOutline(fit$acmap)
  map_ag_names   <- Racmacs::agNames(fit$acmap)
  map_sr_names   <- Racmacs::srNames(fit$acmap)
  
  # Get the plot pars
  pars <- do.call(ablandscape.par, options)
  
  ## Add antigens and sera
  if(show.map.antigens){
    for (n in seq_len(nrow(map_ag_coords))) {
      data3js <- r3js::points3js(data3js,
                                 x          = map_ag_coords[n,1],
                                 y          = map_ag_coords[n,2],
                                 z          = zlim[1]+0.019,
                                 size       = pars$cex.basemap.ags,
                                 col        = map_ag_fill[n],
                                 opacity    = pars$opacity.basemap.ags,
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
                                 z          = zlim[1]+0.019,
                                 pch        = 0,
                                 size       = pars$cex.basemap.sr,
                                 col        = map_sr_outline[n],
                                 opacity    = pars$opacity.basemap.sr,
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
