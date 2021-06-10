
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
predict.ablandscape.fit <- function(
  object,
  coords,
  interval   = "none",
  level      = 0.95,
  negll.fit  = FALSE,
  crop2chull = TRUE,
  control    = list()
  ){
  
  # Take upper and lower bounds for the log titers
  if("ablandscape.delta.fit" %in% class(object)){
    upper_bounds <- object$logtiters.delta.upper
    lower_bounds <- object$logtiters.delta.lower
  } else {
    upper_bounds <- object$logtiters.upper
    lower_bounds <- object$logtiters.lower
  }
  
  # Fitting a vector of titers
  if(is.null(dim(object$titers))){
    
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
    
  } 
  
  # Fitting a matrix of titers
  else {
    
    fit_output <- lapply(seq_len(nrow(object$titers)), function(x){
      
      # Do the fit
      coordfit <- fit_lndscp2coords(
        coords     = coords,
        ag_coords  = object$coords,
        titers     = object$titers[x,],
        max_titers = upper_bounds[x,],
        min_titers = lower_bounds[x,],
        bandwidth  = object$bandwidth,
        degree     = object$degree,
        error.sd   = object$error.sd,
        crop2chull = crop2chull,
        control    = object$control
      )
      
      # Get height and negative log-likelihood of the fits
      list(
        height_estimates = vapply(coordfit, function(x){ x$par[1] }, numeric(1)),
        negll_estimates  = vapply(coordfit, function(x){ x$negll }, numeric(1))
      )
      
    })
    
    individual_height_estimates <- do.call(rbind, lapply(fit_output, function(x){ x$height_estimates }))
    individual_negll_estimates  <- do.call(rbind, lapply(fit_output, function(x){ x$negll_estimates }))
    
    height_estimates <- colMeans(individual_height_estimates)
    negll_estimates  <- colMeans(individual_negll_estimates)
    
    # Calculate confidence intervals if required
    if(interval == "none"){
      
      fitted_values <- height_estimates
      
    } else if(interval == "confidence"){
      
      confidence <- apply(
        individual_height_estimates, 2,
        Hmisc::smean.cl.normal,
        conf.int = level
      )
      
      fitted_values <- cbind(
        height_estimates,
        confidence["Lower",],
        confidence["Upper",]
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
    
  }
  
  
  
  # Return the output
  output
  
}

#' @export
predict_lndscp_grid <- function(
  fit,
  grid_x = NULL,
  grid_y = NULL,
  grid_spacing = 0.5,
  crop2chull   = TRUE,
  padding = 1,
  format = "wide"
){
  
  # Decide on default coordinates for grid lines
  if(is.null(grid_x)) grid_x <- seq(from = min(fit$coords[,1])-padding, to = max(fit$coords[,1])+padding, by = grid_spacing)
  if(is.null(grid_y)) grid_y <- seq(from = min(fit$coords[,2])-padding, to = max(fit$coords[,2])+padding, by = grid_spacing)
  
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
    fit,
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
  if(format == "wide"){
    lndscp_grid
  } else if(format == "long") {
    data.frame(
      x = grid_coords[,1],
      y = grid_coords[,2],
      z = lndscp_heights
    )
  } else {
    stop("'format' must be one of 'wide' or 'long'")
  }
  
}



