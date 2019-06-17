
#' Fitting with antibody landscapes
#' 
#' Fit an antibody landscape to titers and antigen positions.
#'
#' @param titers A named vector of titers against each antigen.
#' @param acmap The antigenic map object containing the antigens titrated.
#' @param bandwidth The bandwidth of the local regression.
#' @param degree The degree of the polynomials to be used, normally 1 or 2.
#' @param control Control parameters: see (\code{\link{ablandscape.control}})
#' @param grid_spacing Grid spacing for surface interpolation
#'
#' @return
#' @export
#'
ablandscape <- function(titers,
                        acmap,
                        bandwidth,
                        degree,
                        grid_spacing = 1,
                        control = list()){
  
  # Match up titers and the map
  ag_matches <- match(names(titers), acmap$ag_names)
  if(sum(is.na(ag_matches)) > 0){
    stop("Not all titer names were present in the antigenic map")
  }
  
  # Make the fit object
  object <- ablandscape.fit(titers    = titers,
                            coords    = acmap$ag_coords[ag_matches,,drop=FALSE],
                            bandwidth = bandwidth,
                            degree    = degree,
                            control   = control)
  
  # Update the class of the object
  class(object) <- c("ablandscape", "ablandscape.fit", "list")
  
  # Build the landscape fit
  object$landscape  <- build_lndscp(object, grid_spacing = grid_spacing)
  object$acmap      <- acmap
  object$ag_indices <- ag_matches
  
  # Return the ablandscape object
  object
  
}





