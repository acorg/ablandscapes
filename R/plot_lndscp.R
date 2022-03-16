
#' View method for the antibody landscape object
#'
#' @param object The antibody landscape object
#' @param xlim The x axis limits
#' @param ylim The y axis limits
#' @param zlim The z axis limits
#' @param crop2chull Should the plot be cropped to the convex hull bounded by the points
#' @param ... Additional plotting parameters
#'
#' @export
#' 
lndscp_r3jsplot <- function(
  fit,
  xlim                        = NULL, 
  ylim                        = NULL, 
  zlim                        = NULL,
  padding                     = 1,
  crop2chull                  = TRUE,
  output                      = "widget",
  show.titers                 = TRUE,
  show.impulses               = show.titers,
  show.surface                = TRUE,
  show.individual.surfaces    = FALSE,
  show.map.antigens           = TRUE,
  show.map.sera               = TRUE,
  show.axis                   = TRUE,
  show.sidegrid               = TRUE,
  aspect.z                    = 1,
  rotation                    = NULL,
  translation                 = NULL,
  zoom                        = NULL,
  toggles                     = TRUE,
  grid_spacing                = 0.5,
  titer_proportion_cutoff     = 0.9,
  options                     = list(),
  options.individual.surfaces = list(
    opacity.surface.grid = 0.4,
    opacity.surface      = 0.2,
    col.surface          = "grey70",
    col.surface.grid     = "grey70"
  )
){
  
  ## Work out plot limits
  if(is.null(fit$acmap)){
    map_coords <- fit$coords
  } else {
    map_coords <- rbind(
      Racmacs::agCoords(fit$acmap),
      Racmacs::srCoords(fit$acmap)
    )
  }
  if(is.null(xlim)) xlim <- calc_map_lims(range(map_coords[,1]), padding = padding)
  if(is.null(ylim)) ylim <- calc_map_lims(range(map_coords[,2]), padding = padding)
  if(is.null(zlim)) zlim <- c(-1, ceiling(max(fit$logtiters, na.rm = TRUE)))
  
  # Setup the plot
  data3js <- lndscp3d_setup(
    xlim     = xlim, 
    ylim     = ylim, 
    zlim     = zlim,
    aspect.z = aspect.z,
    show.axis = show.axis,
    show.sidegrid = show.sidegrid,
    options  = options
  )
  
  # Plot basemap
  if(!is.null(fit$acmap)){
    data3js <- lndscp3d_map(
      data3js = data3js,
      fit     = fit, 
      zlim    = zlim,
      show.map.antigens = show.map.antigens,
      show.map.sera = show.map.sera,
      options = options
    )
  }
  
  # Plot titers
  if(show.titers){
    data3js <- lndscp3d_titers(
      data3js       = data3js,
      object        = fit, 
      zlim          = zlim,
      show.impulses = show.impulses,
      titer_proportion_cutoff = titer_proportion_cutoff,
      options = options
    )
  }

  # Plot surface
  if(show.surface){
    data3js <- lndscp3d_surface(
      data3js      = data3js,
      object       = fit,
      crop2chull   = crop2chull,
      grid_spacing = grid_spacing,
      options = options
    )
  }
  
  if(show.individual.surfaces){
    for (i in seq_len(nrow(fit$titers))) {
      
      individual_fit <- fit
      individual_fit$titers <- fit$titers[i,]
      individual_fit$logtiters <- fit$logtiters[i,]
      individual_fit$logtiters.upper <- individual_fit$logtiters.upper[i,]
      individual_fit$logtiters.lower <- individual_fit$logtiters.lower[i,]
      individual_fit$lessthans <- individual_fit$lessthans[i,]
      individual_fit$morethans <- individual_fit$morethans[i,]
      individual_fit$fitted.values <- NULL
      individual_fit$residuals <- NULL
      individual_fit$residuals.lessthan <- NULL
      individual_fit$residuals.morethan <- NULL
      if (!is.null(individual_fit$cone)) {
        individual_fit$cone$cone_coords  <- individual_fit$cone$cone_coords[i,,drop=F]
        individual_fit$cone$cone_heights <- individual_fit$cone$cone_heights[i]
      }
      
      data3js <- lndscp3d_surface(
        data3js      = data3js,
        object       = individual_fit,
        crop2chull   = crop2chull,
        grid_spacing = grid_spacing,
        options      = options.individual.surfaces,
        toggle       = "Individual surfaces"
      )
      
    }
  }
  
  # Return 3js data or widget
  if(output == "widget"){
    r3js::r3js(
      data3js,
      rotation    = rotation,
      translation = translation,
      zoom        = zoom
    )
  } else {
    data3js
  }
  
}

