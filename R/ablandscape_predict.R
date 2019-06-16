
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
predict.ablandscape <- function(object,
                                coords,
                                interval = c("none", "confidence"),
                                crop2chull = TRUE){
  
  if("ablandscape.delta" %in% class(object)){
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
  
  fitted_values <- vapply(coordfit, function(x){ x$par[1] }, numeric(1))
  fitted_values
  
}


