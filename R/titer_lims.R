
#' Get titer limits
#'
#' Function for getting upper and lower limits of measured titers on the log scale.
#'
#' @param titers A numeric/character vector or matrix of titer measurements.
#' @param fit_opts A list of fitting options (see details).
#'
#' @return Returns a list of length two with values max_titers and min_titers, giving the 
#' numeric vectors of the upper and lower bounds of the titers on the log scale.
#'
#' @details This function assumes that HI measurements were performed in 2-fold dilution 
#' steps and converts them to the log scale using the formula:
#' \deqn{a + b}
#' Hence an HI titer of 20, which would convert to 1 via the transformation above, would be 
#' assumed to have upper and lower limits of 1.5 and 0.5 respectively.
#' 
#' In the case of non-detectable titers, such as <10, the lower bound of the measured value 
#' is taken from the parameter `min_titer_possible`, defaulting to the value found from 
#' a call to `ablandscape.control()`. For a greater than value, i.e. >1280, the 
#' upper bound of the value is taken from the parameter `max_titer_possible`. You can 
#' set different defaults by passing them as named arguments to the list, as shown in the 
#' examples.
#' 
calc_titer_lims <- function(
  titers, 
  control = list()
  ){
  
  # Throw an error if titers is not a matrix or vector
  if(!is.matrix(titers) & !is.vector(titers)){
    stop("Titers must be in the form of a vector or matrix")
  }
  
  # Convert titers to a vector if necessary
  titer_dims <- dim(titers)
  if(is.matrix(titers)){
    titers <- as.vector(titers)
  }
  
  # Find less than and greater than titers and convert them to a numeric form
  lessthan_titers   <- grepl(x = titers, pattern = "<[0-9]")
  lessthaneq_titers <- grepl(x = titers, pattern = "<=[0-9]")
  morethan_titers   <- grepl(x = titers, pattern = ">[0-9]")
  morethaneq_titers <- grepl(x = titers, pattern = ">=[0-9]")
  na_titers       <- grepl(x = titers, pattern = "\\*")
  
  numeric_titers <- titers
  numeric_titers[na_titers] <- NA
  numeric_titers <- as.numeric(gsub("(<|>|<=|>=)","",numeric_titers))
  
  # Convert titers to the log scale
  log_titers <- log2(numeric_titers/10)
  log_titers[lessthan_titers] <- log_titers[lessthan_titers] - 1
  log_titers[morethan_titers] <- log_titers[morethan_titers] + 1
  max_titers <- log_titers + 0.5
  min_titers <- log_titers - 0.5
  
  pars <- do.call(ablandscape.control, control)
  min_titers[lessthan_titers]   <- pars$min.titer.possible
  min_titers[lessthaneq_titers] <- pars$min.titer.possible
  max_titers[morethan_titers]   <- pars$max.titer.possible
  max_titers[morethaneq_titers] <- pars$max.titer.possible
  
  if(!is.null(titer_dims)){
    max_titers <- matrix(data = max_titers, 
                         nrow = titer_dims[1], 
                         ncol = titer_dims[2])
    min_titers <- matrix(data = min_titers, 
                         nrow = titer_dims[1], 
                         ncol = titer_dims[2])
  }
  
  list(max_titers = max_titers,
       min_titers = min_titers)
  
}


#' Get titer differences
#' 
#' Function for getting upper and lower limits of differences between measured titers on the log scale.
#'
#' @param titers1 A numeric/character vector or matrix of pre titer measurements.
#' @param titers2 A numeric/character vector or matrix of post titer measurements.
#' @param control A list of fitting options (see details).
#'
#' @return Returns a list of length two with values max_diffs and min_diffs, giving the numeric 
#' vectors of the upper and lower bounds of the difference in pre and post titers on the log scale.
#'
calc_titer_diffs <- function(
  titers1,
  titers2,
  control = list()
  ) {
  
  # Get limits to titer measurements
  pre_lims  <- calc_titer_lims(titers1, control)
  post_lims <- calc_titer_lims(titers2, control)
  
  # Compute maximum and mininmum differences
  max_diffs <- post_lims$max_titers - pre_lims$min_titers
  min_diffs <- post_lims$min_titers - pre_lims$max_titers
  
  # Return result
  list(max_diffs = max_diffs,
       min_diffs = min_diffs)
  
}




