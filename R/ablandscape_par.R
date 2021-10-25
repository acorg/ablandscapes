
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
  cex.basemap.ags      = cex.basemap,
  cex.basemap.sr       = cex.basemap,
  cex.titer            = cex,
  col.impulse          = "grey20",
  lwd.impulse          = 0.5,
  asp.z                = 1,
  zaxt                 = "linear",
  col.grid             = "grey80",
  col.surface          = "grey20",
  col.surface.grid     = "grey20",
  shininess.surface    = 30,
  opacity.basemap      = 0.4,
  opacity.basemap.ags  = opacity.basemap,
  opacity.basemap.sr   = opacity.basemap,
  opacity.surface      = 0.8,
  opacity.surface.grid = 0.6
) {
  
  list(
    cex                  = cex,
    cex.basemap          = cex.basemap,
    cex.basemap.ags      = cex.basemap.ags,
    cex.basemap.sr       = cex.basemap.sr,
    cex.titer            = cex.titer,
    asp.z                = asp.z,
    zaxt                 = zaxt,
    col.impulse          = col.impulse,
    lwd.impulse          = lwd.impulse,
    col.grid             = col.grid,
    col.surface          = col.surface,
    col.surface.grid     = col.surface.grid,
    opacity.basemap      = opacity.basemap,
    opacity.basemap.ags  = opacity.basemap.ags,
    opacity.basemap.sr   = opacity.basemap.sr,
    opacity.surface      = opacity.surface,
    opacity.surface.grid = opacity.surface.grid,
    shininess.surface    = shininess.surface
  )
  
}


