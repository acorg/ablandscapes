
#' Add a surface to a 3d landscape plot
#'
#' @param data3js The r3js data object
#' @param object The antibody landscape object
#' @param ... Additional plotting parameters
#' 
#' @export
#' 
lndscp3d_surface <- function(data3js,
                             object,
                             crop2chull,
                             ...) {
  
  # Get the plot pars
  pars <- ablandscape.par(...)
  
  # Generate clipping planes
  if(crop2chull){
  
    chullpoints <- chull(object$coords)
    chullpoints <- c(chullpoints, chullpoints[1])
    
    clipping_planes <- lapply(seq_len(length(chullpoints)-1), function(x){

      r3js::clippingPlane3js(coplanarPoints = rbind(
        c(object$coords[chullpoints[x],],   0),
        c(object$coords[chullpoints[x+1],], 0),
        c(object$coords[chullpoints[x+1],], 1)
      ))
  
    })
    
  } else {
  
    clipping_planes <- NULL
    
  }
  
  # Add the surface
  data3js <- r3js::surface3js(data3js,
                              x = object$landscape$x,
                              y = object$landscape$y,
                              z = object$landscape$z,
                              col        = pars$col.surface,
                              toggle     = "Surface",
                              opacity    = pars$opacity.surface,
                              depthWrite = FALSE,
                              xpd        = FALSE,
                              clippingPlanes = clipping_planes,
                              ...)
  
  # Add the grid
  data3js <- r3js::surface3js(data3js,
                              x = object$landscape$x,
                              y = object$landscape$y,
                              z = object$landscape$z,
                              col        = pars$col.surface.grid,
                              wireframe  = TRUE,
                              toggle     = "Surface",
                              opacity    = pars$opacity.surface.grid,
                              depthWrite = FALSE,
                              xpd        = FALSE,
                              clippingPlanes = clipping_planes,
                              ...)
  
  # Return the surface button data
  invisible(data3js)
  
}



