

#' Convert raw titers to log titers
#'
#' @param titers The raw titers to convert.
#'
#' @export
#'
as.logtiter <- function(titers){
  
  if(!is.matrix(titers)){
    titers <- as.matrix(titers)
  }
  mode(titers) <- "character"
  converted_titers <- convert2logCpp(titers)
  converted_titers <- lapply(converted_titers, function(x){
    colnames(x) <- colnames(titers)
    rownames(x) <- rownames(titers)
    x
  })
  converted_titers$log_titers
  
}







