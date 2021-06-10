
calc_map_lims <- function(lims, padding = 1, round_even = TRUE){
  lims[1] <- lims[1] - padding
  lims[2] <- lims[2] + padding
  if(round_even){
    d  <- diff(lims)
    dd <- ceiling(d) - d
    lims[1] <- lims[1] - dd/2
    lims[2] <- lims[2] + dd/2
  }
  lims
}

in_range <- function(x, coords){
  x > min(coords) && x <= max(coords)
}

cross_point_x <- function(l, x){
  l <- l[order(l[,2]),]
  gradient <- diff(l[,2]) / diff(l[,1])
  diffx <- x - l[1,1]
  l[1,2] + gradient*diffx
}

cross_point_y <- function(l, y){
  l <- l[order(l[,1]),]
  gradient <- diff(l[,1]) / diff(l[,2])
  diffy <- y - l[1,2]
  l[1,1] + gradient*diffy
}

#' Plot an ablandscape contour plot
#' 
#' Plot a contour plot from an antibody landscape fit using `ggplot2` and [ggplot2::geom_contour_filled()].
#' 
#' @param fit The antibody landscape fit object, created using [ablandscape.fit()]
#' @param xlim x limits for the plot
#' @param ylim y limits for the plot
#' @param crop2chull Crop the plot to the convex hull defined by the antigens titrated
#' @param padding Padding to apply to the fit
#' @param grid_x x grid points at which the ablandscapes fit should be calculated (passed to [predict_lndscp_grid()])
#' @param grid_y y grid points at which the ablandscapes fit should be calculated (passed to [predict_lndscp_grid()])
#' @param grid_spacing grid spacing for the ablandscapes fit (passed to [predict_lndscp_grid()])
#' @param show.grid Show grid lines behind contours
#' @param grid.color Grid line colors when shown
#' @param grid.linewidth Grid line width when shown
#' @param ... Additional arguments passed to [ggplot2::geom_contour_filled()]
#' 
#' @return Returns the ggplot object for the contour plot
#' 
#' @export
lndscp_contour_ggplot <- function(
  fit,
  xlim = NULL,
  ylim = NULL,
  crop2chull = TRUE,
  padding = 1,
  grid_x = NULL,
  grid_y = NULL,
  grid_spacing = 0.5,
  show.grid = TRUE,
  grid.color = "grey90",
  grid.linewidth = 0.5,
  ...
){
  
  # Check for required packages
  if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("Package \"geometry\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Fit the grid
  lndscp_grid <- predict_lndscp_grid(
    fit, 
    crop2chull = FALSE,
    format = "long",
    grid_x = grid_x,
    grid_y = grid_y,
    grid_spacing = grid_spacing,
    padding = padding
  )
  
  # Define the fit coords
  fit_coords <- fit$coords
  
  # Set default plot limits
  if(crop2chull){
    if(is.null(xlim)) xlim <- calc_map_lims(range(fit_coords[,1]), padding = padding)
    if(is.null(ylim)) ylim <- calc_map_lims(range(fit_coords[,2]), padding = padding)
  } else {
    if(is.null(xlim)) xlim <- range(lndscp_grid$x)
    if(is.null(ylim)) ylim <- range(lndscp_grid$y)
  }
  
  # Do the basic contour plot
  ggplot2::ggplot(
    lndscp_grid
  ) + 
    ggplot2::geom_contour_filled(
      ggplot2::aes(
        x = x,
        y = y,
        z = z
      ),
      ...
    ) + 
    ggplot2::coord_fixed(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) + 
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA
      )
    ) -> gp
  
  # Mask regions outside the convex hull of points titrated
  if(crop2chull){
    
    # Get fit coords and perform some delauney triangulation
    # to work out masking polygons
    coords <- fit_coords[chull(fit_coords),]
    coords <- rbind(
      c(xlim[1], ylim[1]),
      c(xlim[1], ylim[2]),
      c(xlim[2], ylim[2]),
      c(xlim[2], ylim[1]),
      coords
    )
    
    dcoords <- geometry::delaunayn(coords)
    dcoords <- dcoords[apply(dcoords, 1, function(x){ sum(x %in% 1:4) > 0 }),]
    
    # Add the masking polygons
    for(x in seq_len(nrow(dcoords))) {
      gp <- gp + ggplot2::annotate(
        geom = "polygon",
        x = coords[,1][dcoords[x,]],
        y = coords[,2][dcoords[x,]],
        fill = "#ffffff",
        color = "#ffffff"
      )
    }
    
    # Add in grid lines that stop at the masking edge
    if(show.grid){
      
      chull_indices <- chull(fit_coords)
      outlines <- lapply(seq_along(chull_indices), function(x){
        if(x == length(chull_indices)) x2 <- 1
        else                           x2 <- x + 1
        fit_coords[chull_indices[c(x,x2)],]
      })
      
      # Draw horizontal lines up until edge of masked contour
      for(x in seq(from = ylim[1], to = ylim[2])){
        outline_crossed <- vapply(outlines, function(l){
          in_range(x, l[,2])
        }, logical(1))
        if(sum(outline_crossed) == 0){
          gp <- gp + ggplot2::geom_hline(yintercept = x, color = grid.color, size = grid.linewidth)
        } else {
          cross_points <- vapply(
            outlines[outline_crossed],
            cross_point_y,
            numeric(1),
            y = x
          )
          cross_points <- sort(cross_points)
          gp <- gp + ggplot2::annotate("line", x = c(xlim[1], cross_points[1]), y = rep(x,2), color = grid.color, size = grid.linewidth)
          gp <- gp + ggplot2::annotate("line", x = c(cross_points[2], xlim[2]), y = rep(x,2), color = grid.color, size = grid.linewidth)
        }
      }
      
      # Draw vertical lines up until edge of masked contour
      for(x in seq(from = xlim[1], to = xlim[2])){
        outline_crossed <- vapply(outlines, function(l){
          in_range(x, l[,1])
        }, logical(1))
        if(sum(outline_crossed) == 0){
          gp <- gp + ggplot2::geom_vline(xintercept = x, color = grid.color, size = grid.linewidth)
        } else {
          cross_points <- vapply(
            outlines[outline_crossed],
            cross_point_x,
            numeric(1),
            x = x
          )
          cross_points <- sort(cross_points)
          gp <- gp + ggplot2::annotate("line", x = rep(x,2), y = c(ylim[1], cross_points[1]), color = grid.color, size = grid.linewidth)
          gp <- gp + ggplot2::annotate("line", x = rep(x,2), y = c(cross_points[2], ylim[2]), color = grid.color, size = grid.linewidth)
        }
      }
      
    }
    
  }
  
  # Return the plot
  gp
  
}


