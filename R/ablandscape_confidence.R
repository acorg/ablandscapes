

#' Get landscape confidence interval
#'
#' Function for calculating confidence intervals.
#'
fit_ci2coords <- function(coords,
                          heights,
                          negll,
                          ag_coords,
                          max_titers,
                          min_titers,
                          bandwidth,
                          degree,
                          crop2chull,
                          control){
  
  
  for(x in 1){
    result <- optim(par        = heights[x],
                    fn         = target_height,
                    method     = "L-BFGS-B",
                    coords     = coords[x,],
                    negll      = negll[x],
                    ag_coords  = ag_coords,
                    max_titers = max_titers,
                    min_titers = min_titers,
                    bandwidth  = bandwidth,
                    degree     = degree,
                    crop2chull = crop2chull,
                    control2   = control)
    
  }
  
}


target_height <- function(height, coords, negll, control2, ...){
  
  x <- abs(negll - fit_height2point(lndscp_height = height, 
                                    point_coords  = coords,
                                    control       = control2,
                                    ...)$negll*-2 - qchisq(0.05, 1))
  print(x)
  
}



#' Fit landscape to single coordinate point
#'
fit_height2point <- function(point_coords,
                             lndscp_height,
                             ag_coords,
                             max_titers,
                             min_titers,
                             bandwidth,
                             degree,
                             crop2chull = TRUE,
                             control    = list()) {
  
  # Return NA if outside convex hull
  if(crop2chull) {
    withinChull <- lndscp_checkChull(point_coords, ag_coords)
    if(!withinChull) {
      return(
        list(par   = c(NA, NA, NA),
             negll = NA)
      )
    }
  }
  
  # Get model variables to send to optimiser
  pars <- do.call(ablandscape.control, control)
  
  # Make ag coords relative to the point coord
  ag_coords <- ag_coords - matrix(point_coords, 
                                  nrow = nrow(ag_coords), 
                                  ncol = ncol(ag_coords),
                                  byrow = TRUE)
  
  # Make coords polynomial
  ag_coords_poly <- polycoords(ag_coords, degree)
  
  # Perform the fit
  result <- nlminb(start         = rep(0, ncol(ag_coords_poly)),
                   objective     = get_negll_hi_height,
                   upper         = rep(pars$max.slope, ncol(ag_coords)),
                   lower         = rep(-pars$max.slope, ncol(ag_coords)),
                   ag_coords     = ag_coords_poly,
                   ag_weights    = agweights(ag_coords, bandwidth),
                   max_titres    = matrix(max_titers, nrow = 1),
                   min_titres    = matrix(min_titers, nrow = 1),
                   error_sd      = pars$error.sd,
                   lndscp_height = lndscp_height)
  
  # Return the best height fit
  list(par   = result$par,
       negll = result$objective)
  
}

target_negll <- function(start,
                         pars){
  
  
  
}















