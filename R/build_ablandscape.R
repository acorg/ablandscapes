
#' Build an antibody landscape
#'
#' This is a function for building and antibody landscape from a given set of antigenic coordinates and titers.
#'
#' @param object 
#' @param grid_spacing 
#'
#' @return Returns a list, with a matrix of x and y coordinates for each point
#'   to be plotted on the landscape, and a matrix of z coordinates, which are
#'   the landscape heights at each of these points.
#'
build_lndscp <- function(object,
                         grid_spacing = 1) {
  
  # Grid the coordinates
  coordgrid <- grid_coords(coords       = object$coords,
                           grid_spacing = grid_spacing)
  
  # Fit the z coordinates
  coordgrid$z[] <- predict(object,
                           coords = cbind(
                             as.vector(coordgrid$x),
                             as.vector(coordgrid$y)
                           ),
                           crop2chull = FALSE)
  
  # Return the coordinate grid
  coordgrid
  
}


# Generate an evenly space grid over the coordinates
grid_coords <- function(coords, grid_spacing = 1){
  
  xlines <- seq(from = min(coords[,1]) - grid_spacing/2,
                to   = max(coords[,1]) + grid_spacing/2,
                by   = grid_spacing)
  
  ylines <- seq(from = min(coords[,2]) - grid_spacing/2,
                to   = max(coords[,2]) + grid_spacing/2,
                by   = grid_spacing)
  
  list(
    x = matrix(data  = rep(xlines, times = length(ylines)), 
               ncol  = length(xlines), 
               nrow  = length(ylines), 
               byrow = T),
    y = matrix(data  = rep(ylines, each  = length(xlines)), 
               ncol  = length(xlines), 
               nrow  = length(ylines), 
               byrow = T),
    z = matrix(ncol  = length(xlines), 
               nrow  = length(ylines), 
               byrow = T)
  )
  
}

