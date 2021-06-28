
#' Get defaults for antibody landscape fitting
#'
#' This function returns all the options that can be set for fitting of antibody
#' landscapes and the default values for them.
#'
#' @return Returns a list of default fitting options for making antibody
#'   landscapes.
#' 
#' @export
#' 
ablandscape.control <- function(
  max.titer.possible = 10,
  min.titer.possible = -4,
  max.slope          = Inf,
  error.sd           = 1.1,
  confint.stepsize   = 0.1,
  model.fn           = NULL
) {
  
  
  list(
    max.titer.possible = max.titer.possible,
    min.titer.possible = min.titer.possible,
    max.slope          = max.slope,
    error.sd           = error.sd,
    confint.stepsize   = confint.stepsize,
    model.fn           = model.fn
  )
  
}

