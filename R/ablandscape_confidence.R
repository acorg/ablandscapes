
#' Get landscape confidence interval
#'
#' Function for calculating confidence intervals.
#'
fit_ci2coords <- function(
  coords,
  heights,
  negll,
  ag_coords,
  max_titers,
  min_titers,
  bandwidth,
  degree,
  error.sd,
  crop2chull,
  control,
  prediction_control,
  level
  ) {
  
  # Get model variables for confidence intervals
  pars            <- do.call(ablandscape.control, control)
  prediction_pars <- do.call(ablandscape.control, prediction_control)
  
  stepsize <- prediction_pars$confint.stepsize
  
  # Setup to store upper confidence intervals and lower confidence intervals
  upper_ci <- rep(NA, length(heights))
  lower_ci <- rep(NA, length(heights))
  
  # Work out upper confidence intervals for each coordinate
  for (i in seq_along(heights)) {
    
    # Calculate the target negll for the upper and lower confidence interval
    target_negll <- negll[i] + qchisq(level, 1)/2
    
    # Start from the minimum negll
    upper_negll  <- negll[i]
    upper_height <- heights[i]
    
    # Increase the landscape height until the target negll is reached
    while (upper_negll < target_negll) {
      upper_negll_last  <- upper_negll
      upper_height_last <- upper_height
      upper_height <- upper_height + stepsize
      upper_negll  <- fit_height2point(
        lndscp_height = upper_height, 
        point_coords  = coords[i,],
        control       = control,
        ag_coords     = ag_coords,
        max_titers    = max_titers,
        min_titers    = min_titers,
        bandwidth     = bandwidth,
        degree        = degree,
        error.sd      = error.sd,
        crop2chull    = crop2chull
      )$negll
    }
    
    # Work out where upperci should roughly fall between stepsizes
    upper_ci[i] <- upper_height_last +
                     (upper_height - upper_height_last) *
                     ((target_negll - upper_negll_last) / (upper_negll - upper_negll_last))
    
    # Start from the minimum negll
    lower_negll  <- negll[i]
    lower_height <- heights[i]
    
    # Decrease the landscape height until the target negll is reached
    while (lower_negll < target_negll) {
      lower_negll_last  <- lower_negll
      lower_height_last <- lower_height
      lower_height <- lower_height - stepsize
      lower_negll  <- fit_height2point(
        lndscp_height = lower_height, 
        point_coords  = coords[i,],
        control       = control,
        ag_coords     = ag_coords,
        max_titers    = max_titers,
        min_titers    = min_titers,
        bandwidth     = bandwidth,
        degree        = degree,
        error.sd      = error.sd,
        crop2chull    = crop2chull
      )$negll
    }
    
    # Work out where upperci should roughly fall between stepsizes
    lower_ci[i] <- lower_height_last +
      (lower_height - lower_height_last) *
      ((target_negll - lower_negll_last) / (lower_negll - lower_negll_last))
    
  }
  
  # Return the confidence intervals
  cbind(lower_ci, upper_ci)
  
}




#' Fit landscape to single coordinate point
#'
fit_height2point <- function(
  point_coords,
  lndscp_height,
  ag_coords,
  max_titers,
  min_titers,
  bandwidth,
  degree,
  error.sd,
  crop2chull = TRUE,
  control    = list()
  ) {
  
  # Return NA if outside convex hull
  if (crop2chull) {
    withinChull <- lndscp_checkChull(point_coords, ag_coords)
    if (!withinChull) {
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
  result <- nlminb(
    start         = rep(0, ncol(ag_coords_poly)),
    objective     = negll_lndscp_height,
    upper         = rep(pars$max.slope, ncol(ag_coords)),
    lower         = rep(-pars$max.slope, ncol(ag_coords)),
    ag_coords     = ag_coords_poly,
    ag_weights    = tricubic.weights(ag_coords, bandwidth),
    max_titers    = matrix(max_titers, nrow = 1),
    min_titers    = matrix(min_titers, nrow = 1),
    error_sd      = error.sd,
    lndscp_height = lndscp_height
  )
  
  # Return the best height fit
  list(
    par   = result$par,
    negll = result$objective
  )
  
}
