
lndscp3d_setup <- function(
  xlim,
  ylim,
  zlim,
  aspect.z,
  show.sidegrid = TRUE,
  show.axis = TRUE,
  options
){
  
  # Get the plot pars
  pars <- do.call(ablandscape.par, options)

  # Open a new r3js plot
  data3js <- r3js::plot3js.new()
  
  ## Set plot window
  data3js <- r3js::plot3js.window(
    data3js,
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    aspect = c(1, 1, aspect.z)
  )
  
  ## Set box
  if (show.sidegrid) {
    data3js <- r3js::box3js(
      data3js,
      col   = "grey80"
    )
  }
  
  ## Add a side grid
  if (show.sidegrid) {
    data3js <- r3js::grid3js(
      data3js,
      axes = "z",
      lwd = pars$sidegrid.lwd,
      col = pars$sidegrid.col,
      at = pars$sidegrid.at
    )
  }
  
  ## Display axes
  if (show.axis) {
    axis_at <- zlim[1]:zlim[2]
    if (pars$zaxt == "linear") {
      axis_labels <- 2^axis_at*10
      axis_labels[axis_labels == 5] <- "<10"
    }
    if (pars$zaxt == "log") {
      axis_labels <- axis_at
      axis_labels[axis_labels == -1] <- "nd"
    }
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      cornerside = "f"
    )
  }
  
  ## Add a base grid
  for(x in seq(from = xlim[1], to = xlim[2])){
    data3js <- r3js::lines3js(
      data3js,
      x = c(x, x),
      y = range(ylim),
      z = c(zlim[1], zlim[1]),
      col = pars$col.grid,
      lwd = pars$lwd.grid,
      xpd = TRUE
    )
  }
  for(y in seq(from = ylim[1], to = ylim[2])){
    data3js <- r3js::lines3js(
      data3js,
      x = range(xlim),
      y = c(y, y),
      z = c(zlim[1], zlim[1]),
      col = pars$col.grid,
      lwd = pars$lwd.grid,
      xpd = TRUE
    )
  }
  
  # Return the data
  data3js
  
}

