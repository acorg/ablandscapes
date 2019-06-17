
#' Get defaults for antibody landscape plotting
#'
#' This function returns all the options that can be set for plotting of antibody
#' landscapes and the default values for them.
#'
#' @return Returns a list of default plotting options for making antibody landscapes.
#' 
#' @export
#' 
ablandscape.par <- function(
  cex                  = 1,
  cex.basemap          = cex,
  cex.titer            = cex*2,
  asp.z                = 1,
  zaxt                 = "linear",
  col.grid             = "grey80",
  col.surface          = "grey20",
  col.surface.grid     = "grey20",
  opacity.basemap      = 0.8,
  opacity.surface      = 0.8,
  opacity.surface.grid = 0.6
) {
  
  list(
    cex                  = cex,
    cex.basemap          = cex.basemap,
    cex.titer            = cex.titer,
    asp.z                = asp.z,
    zaxt                 = zaxt,
    col.grid             = col.grid,
    col.surface          = col.surface,
    col.surface.grid     = col.surface.grid,
    opacity.basemap      = opacity.basemap,
    opacity.surface      = opacity.surface,
    opacity.surface.grid = opacity.surface.grid
  )
  
}


