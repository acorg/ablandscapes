
#' View objects interactively
#'
#' Generic function for viewing objects interactively.
#'
#' @param x The object
#' @param ... Objects to pass to the method
#'
#' @export
#'
view <- function (x, ...) {
  UseMethod("view", x)
}

#' View method for the antibody landscape object
#'
#' @param object The antibody landscape object
#' @param xlim The x axis limits
#' @param ylim The y axis limits
#' @param zlim The z axis limits
#' @param crop2chull Should the plot be cropped to the convex hull bounded by the points
#' @param ... Additional plotting parameters
#'
#' @export
#' 
view.ablandscape <- function(object,
                             xlim, ylim, zlim,
                             crop2chull        = TRUE,
                             output            = "widget",
                             show.titers       = TRUE,
                             show.impulses     = show.titers,
                             show.surface      = TRUE,
                             show.map.antigens = TRUE,
                             show.map.sera     = TRUE,
                             aspect.z          = 1,
                             rotation          = NULL,
                             translation       = NULL,
                             zoom              = NULL,
                             toggles           = TRUE,
                             ...){
  
  # Get the map data
  map_ag_coords <- Racmacs::agCoords(object$acmap)
  map_sr_coords <- Racmacs::srCoords(object$acmap)
  
  # Get all base coords
  base_coords <- rbind(
    map_ag_coords, 
    map_sr_coords
  )
  
  ## Work out plot limits
  if(missing(xlim)){
    xlim <- c(floor(min(base_coords[,1]))-1,
              ceiling(max(base_coords[,1]))+1)
  }
  if(missing(ylim)){
    ylim <- c(floor(min(base_coords[,2]))-1,
              ceiling(max(base_coords[,2]))+1)
  }
  if(missing(zlim)){
    zlim <- c(-1, ceiling(max(object$logtiters, na.rm = TRUE)))
  }
  
  # Plot basemap
  data3js <- lndscp3d_map(object   = object, 
                          xlim     = xlim, 
                          ylim     = ylim, 
                          zlim     = zlim,
                          aspect.z = aspect.z,
                          show.map.antigens = show.map.antigens,
                          show.map.sera = show.map.sera,
                          toggle = toggles,
                          ...)
  
  # Plot titers
  if(show.titers){
    data3js <- lndscp3d_titres(data3js       = data3js,
                               object        = object, 
                               zlim          = zlim,
                               show.impulses = show.impulses,
                               toggle        = toggles,
                               ...)
  }

  # Plot surface
  if(show.surface){
    data3js <- lndscp3d_surface(data3js    = data3js,
                                object     = object,
                                crop2chull = crop2chull,
                                toggle     = toggles,
                                ...)
  }
  
  # Add an axis
  tickpoints <- seq(from = min(zlim), to = max(zlim), by = 1)
  ticklabels <- 2^tickpoints*10
  ticklabels[ticklabels == 5] <- "<10"
  
  data3js <- r3js::axis3js(
    data3js,
    side = "z",
    cornerside = "f",
    at     = tickpoints,
    labels = ticklabels,
    lwd = 1
  )
  
  # Return 3js data or widget
  if(output == "widget"){
    r3js::r3js(
      data3js,
      rotation    = rotation,
      translation = translation,
      zoom        = zoom
    )
  } else {
    data3js
  }
  
}

