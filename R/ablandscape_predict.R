
#' Predict method for antibody landscapes
#' 
#' Produce predicted values based on a antibody landscapes fit object.
#'
#' @param object The antibody landscapes fit object
#' @param coords Antigenic coordinates for which to predict landscapes heights
#' @param interval Type of interval calculation, (one of "none" or "confidence")
#'
#' @return Produces a vector of predictions or a matrix of predictions and
#'   bounds with column names fit, lwr, and upr if interval is set.
#' @export
#'
predict.ablandscape.fit <- function(object,
                                    coords,
                                    interval   = "none",
                                    level      = 0.95,
                                    negll.fit  = FALSE,
                                    crop2chull = TRUE,
                                    control    = list()){
  
  # Take upper and lower bounds for the log titers
  if("ablandscape.delta.fit" %in% class(object)){
    upper_bounds <- object$logtiters.delta.upper
    lower_bounds <- object$logtiters.delta.lower
  } else {
    upper_bounds <- object$logtiters.upper
    lower_bounds <- object$logtiters.lower
  }
  
  # Fit the new coordinates supplied
  coordfit <- fit_lndscp2coords(
    coords     = coords,
    ag_coords  = object$coords,
    titers     = object$titers,
    max_titers = upper_bounds,
    min_titers = lower_bounds,
    bandwidth  = object$bandwidth,
    degree     = object$degree,
    error.sd   = object$error.sd,
    crop2chull = crop2chull,
    control    = object$control
  )
  
  # Get height and negative log-likelihood of the fits
  height_estimates <- vapply(coordfit, function(x){ x$par[1] }, numeric(1))
  negll_estimates  <- vapply(coordfit, function(x){ x$negll }, numeric(1))
  
  # Calculate confidence intervals if required
  if(interval == "none"){
    
    fitted_values <- height_estimates
    
  } else if(interval == "confidence"){
    
    confidence <- fit_ci2coords(
      coords     = coords,
      heights    = height_estimates,
      negll      = negll_estimates,
      ag_coords  = object$coords,
      max_titers = upper_bounds,
      min_titers = lower_bounds,
      bandwidth  = object$bandwidth,
      degree     = object$degree,
      error.sd   = object$error.sd,
      crop2chull = crop2chull,
      control    = object$control,
      prediction_control = control,
      level      = level
    )
    
    fitted_values <- cbind(
      height_estimates,
      confidence
    )
    
    colnames(fitted_values) <- c("fit", "lwr", "upr")
    
  } else {
    
    stop("Interval must be one of 'none' or 'prediction'")
    
  }
  
  # Return negll of the fit if required
  if(negll.fit){
    output <- list(
      fit       = fitted_values,
      negll.fit = negll_estimates
    )
  } else {
    output <- fitted_values
  }
  
  # Return the output
  output
  
}

#' @export
predict_lndscp_grid <- function(
  object,
  grid_x = NULL,
  grid_y = NULL,
  grid_spacing = 0.5,
  crop2chull   = TRUE
){
  
  # Decide on default coordinates for grid lines
  if(is.null(grid_x)) grid_x <- seq(from = min(object$coords[,1])-1, to = max(object$coords[,1])+1, by = grid_spacing)
  if(is.null(grid_y)) grid_y <- seq(from = min(object$coords[,2])-1, to = max(object$coords[,2])+1, by = grid_spacing)
  
  # Make the grid
  lndscp_grid <- list(
    x = grid_x,
    y = grid_y
  )
  
  # Predict the landscape heights along the grid
  grid_coords    <- expand.grid(
    lndscp_grid$x, 
    lndscp_grid$y
  )
  
  lndscp_heights <- predict(
    object,
    coords     = grid_coords,
    crop2chull = crop2chull
  )
  
  # Convert the predictions back to a matrix matching x and y
  lndscp_grid$z <- matrix(
    lndscp_heights,
    nrow = length(lndscp_grid$x),
    ncol = length(lndscp_grid$y)
  )
  
  # Return the grid object
  lndscp_grid
  
}



