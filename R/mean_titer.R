
#' Calculate the geometric mean titer
#' 
#' This is the basic function for calculating the geometric mean titer of a 
#' set of titers based on maximum likelihood.
#'
#' @param titers A character vector or matrix of HI titers.
#' @param titer_sd The expected standard deviation of the titers around the mean, see details.
#' @param control Control parameters: see (\code{\link{ablandscape.control}})
#'
#' @return Returns the maximum likelihood geometric mean titer.
#' @export
#'
mean.titer <- function(titers,
                       titer_sd = NULL,
                       control = list()){
  
  # Remove NA titers
  if(sum(!is.na(titers)) == 0){ return(NA) }
  titers <- titers[!is.na(titers)]
  
  # Find the titer limits
  titer_lims <- calc_titer_lims(titers  = titers,
                                control = control)
  
  # Calculate the maximum likelihood mean
  mean_fromlims(upper_lims = titer_lims$max_titers,
                lower_lims = titer_lims$min_titers,
                measurement_sd = titer_sd)
  
}


#' Function for calculating the mean foldchange from a group of titers
#'
#' @param titers1 A character vector of starting titers
#' @param titers2 A character vector of end titers
#' @param titer_sd The expected standard deviation of the fold-change responses around the mean, see details.
#' @param control Control parameters: see (\code{\link{ablandscape.control}})
#'
#' @return Returns the maximum likelihood mean fold-change response.
#' @export
#'
mean.foldchange <- function(titers1,
                            titers2,
                            titer_sd = NULL,
                            control = list()){
  
  # Remove NA titers
  na_titers <- is.na(titers1) | is.na(titers2)
  if(sum(!na_titers) == 0){ return(NA) }
  
  titers1  <- titers1[!na_titers]
  titers2 <- titers2[!na_titers]
  
  # Find the titer differences
  titer_diffs <- calc_titer_diffs(titers1  = titers1,
                                  titers2 = titers2,
                                  control     = control)
  
  # Calculate the maximum likelihood mean
  mean_fromlims(upper_lims = titer_diffs$max_diffs,
                lower_lims = titer_diffs$min_diffs,
                measurement_sd = titer_sd,
                start_mean = mean(titer_diffs$max_diffs - 1))
  
}


#' Calculate maximum likelihood mean from measurement limits
#'
#' This is the underlying function called by functions such as mean_titer and 
#' mean_foldchange.
#'
#' @param upper_lims The upper measurement limits
#' @param lower_lims The lower measurement limits
#' @param measurement_sd The expected standard deviation of the measurement accounting for error and variation
#'
#' @return Returns the maximum likelihood mean.
#'
mean_fromlims <- function(upper_lims,
                          lower_lims,
                          measurement_sd = NULL,
                          start_mean = NULL){
  
  # Calculate a starting mean
  if(is.null(start_mean)){
    start_mean <- mean(upper_lims-0.5)
  }
  
  if(is.null(measurement_sd)){
    
    # Set starting parameters
    start_pars   <- c(start_mean, sd(upper_lims))
    
    # Run optimisation
    optim_result <- stats::optim(par = start_pars, 
                                 fn  = calc_mean_titer_negll_without_sd,
                                 max_titers = upper_lims,
                                 min_titers = lower_lims,
                                 method = "L-BFGS-B",
                                 lower = c(min(lower_lims), 0),
                                 upper = c(max(upper_lims), Inf))
    
    # Return the estimate
    estimate <- list(gmt = optim_result$par[1],
                     sd  = optim_result$par[2])
    
  } else {
    
    # Set starting parameters
    start_pars   <- start_mean
    
    # Run optimisation
    optim_result <- stats::optim(par = start_pars, 
                                 fn  = calc_mean_titer_negll,
                                 max_titers = upper_lims,
                                 min_titers = lower_lims,
                                 titer_sd = measurement_sd,
                                 method = "Brent",
                                 lower = min(lower_lims),
                                 upper = max(upper_lims))
    
    # Return the estimate
    estimate <- list(gmt = optim_result$par[1])
    
  }
  
  
  # Return the optim result
  list(estimate  = estimate,
       arguments = list(titer_sd = measurement_sd))
  
}










