
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
                             crop2chull = TRUE,
                             ...){
  
  # Get all base coords
  base_coords <- rbind(
    object$acmap$ag_coords, 
    object$acmap$sr_coords
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
  data3js <- lndscp3d_map(object = object, 
                          xlim   = xlim, 
                          ylim   = ylim, 
                          zlim   = zlim,
                          ...)
  
  # Plot titers
  data3js <- lndscp3d_titres(data3js = data3js,
                             object  = object, 
                             zlim    = zlim, 
                             ...)

  # # Plot surface
  data3js <- lndscp3d_surface(data3js    = data3js,
                              object     = object,
                              crop2chull = crop2chull,
                              ...)
  
  # Return widget
  r3js::r3js(data3js)
  
}

