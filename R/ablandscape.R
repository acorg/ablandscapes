
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
ablandscape.fit_to_acmap <- function(
  titers,
  acmap,
  bandwidth,
  degree,
  error.sd,
  control = list()
  ){
  
  # Check for required packages
  if (!requireNamespace("Racmacs", quietly = TRUE)) {
    stop("Package \"Racmacs\" needed for this function to work. Please install it (see https://acorg.github.io/Racmacs for details).",
         call. = FALSE)
  }
  
  # Get map details
  map_ag_names  <- Racmacs::agNames(acmap)
  map_ag_coords <- Racmacs::agCoords(acmap)
  
  # Match up titers and the map
  ag_matches <- match(names(titers), map_ag_names)
  
  # Warn and ignore titers not found in the map
  ag_mismatches <- is.na(ag_matches)
  if(sum(!ag_mismatches) == 0){
    stop("No matching antigens were found in the map.", call. = FALSE)
  }
  
  if(sum(ag_mismatches) > 0){
    warning(
      "The following strains were not found in the map and have been ignored:\n\n'",
      paste(names(titers)[ag_mismatches], collapse = "'\n'"), "'",
      call. = FALSE
    )
    titers     <- titers[!ag_mismatches]
    ag_matches <- ag_matches[!ag_mismatches]
  }
  
  # Make the fit object
  object <- ablandscape.fit(
    titers    = titers,
    coords    = map_ag_coords[ag_matches,,drop=FALSE],
    bandwidth = bandwidth,
    degree    = degree,
    error.sd  = error.sd,
    control   = control
  )
  
  # Update the class of the object
  class(object) <- c("ablandscape", "ablandscape.fit", "list")
  
  # Build the landscape fit
  object$acmap      <- acmap
  object$ag_indices <- ag_matches
  
  # Return the ablandscape object
  object
  
}





