

#' Convert raw titers to log titers
#'
#' @param titers The raw titers to convert.
#'
#' @export
#'
as.logtiter <- function(titers){
  
  # Check whether a data frame was supplied
  df <- is.data.frame(titers)
  
  # Convert to a matrix
  if(!is.matrix(titers)){
    titers <- as.matrix(titers)
  }
  
  # Make sure the format is character not numeric
  mode(titers) <- "character"
  
  # Perform the conversion
  logtiters <- convert2logCpp(titers)
  
  # Matchup column and row names
  colnames(logtiters) <- colnames(titers)
  rownames(logtiters) <- rownames(titers)
  
  # Convert back to a data frame if necessary
  if(df){
    logtiters <- as.data.frame(logtiters)
  }
  
  # Return the log titers
  logtiters
  
}





