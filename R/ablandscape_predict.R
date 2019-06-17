
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
                                    crop2chull = TRUE){
  
  if("ablandscape.delta.fit" %in% class(object)){
    upper_bounds <- object$logtiters.delta.upper
    lower_bounds <- object$logtiters.delta.lower
  } else {
    upper_bounds <- object$logtiters.upper
    lower_bounds <- object$logtiters.lower
  }
  
  coordfit <- fit_lndscp2coords(
    coords     = coords,
    ag_coords  = object$coords,
    max_titers = upper_bounds,
    min_titers = lower_bounds,
    bandwidth  = object$bandwidth,
    degree     = object$degree,
    crop2chull = crop2chull,
    control    = object$control
  )
  
  height_estimates <- vapply(coordfit, function(x){ x$par[1] }, numeric(1))
  negll_estimates  <- vapply(coordfit, function(x){ x$negll }, numeric(1))
  
  if(interval == "confidence"){
    
    stop("Fitting of confidence intervals is not yet supported, but should be soon. :)")
    # confint <- fit_ci2coords(
    #   coords     = coords,
    #   heights    = height_estimates,
    #   negll      = negll_estimates,
    #   ag_coords  = object$coords,
    #   max_titers = upper_bounds,
    #   min_titers = lower_bounds,
    #   bandwidth  = object$bandwidth,
    #   degree     = object$degree,
    #   crop2chull = crop2chull,
    #   control    = object$control
    # )
    
  }
  
  fitted_values <- height_estimates
  fitted_values
  
}


