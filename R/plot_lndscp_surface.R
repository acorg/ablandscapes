
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
                             crop2chull = TRUE,
                             toggle = "Surface",
                             grid_spacing = 1,
                             padding = 0,
                             options = list()) {
  
  # Get the plot pars
  pars <- do.call(ablandscape.par, options)
  
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
  
  # Predict the grid surface if not done already
  if(is.null(object$grid)){
    object$grid <- predict_lndscp_grid(
      fit = object,
      crop2chull = FALSE,
      padding = padding,
      grid_spacing = grid_spacing
    )
  }
  
  # Add the surface
  data3js <- r3js::surface3js(
    data3js,
    x = object$grid$x,
    y = object$grid$y,
    z = object$grid$z,
    col        = pars$col.surface,
    toggle     = toggle,
    opacity    = pars$opacity.surface,
    shininess  = pars$shininess.surface,
    xpd        = TRUE,
    clippingPlanes = clipping_planes,
    doubleSide = TRUE
  )
  
  # Add the grid
  data3js <- r3js::surface3js(
    data3js,
    x = object$grid$x,
    y = object$grid$y,
    z = object$grid$z,
    col        = pars$col.surface.grid,
    wireframe  = TRUE,
    toggle     = toggle,
    opacity    = pars$opacity.surface.grid,
    xpd        = TRUE,
    clippingPlanes = clipping_planes
  )
  
  # Return the surface button data
  invisible(data3js)
  
}



