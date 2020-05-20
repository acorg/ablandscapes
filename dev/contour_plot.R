
pdf("~/Desktop/test.pdf", 6, 20)
plot.new()
plot.window(
  xlim = xlim,
  ylim = ylim,
  xaxs = "i",
  yaxs = "i",
  asp = 1
)

lndscp_levels <- seq(
  from = floor(min(lndscp_grid$z, na.rm = T)), 
  to   = ceiling(max(lndscp_grid$z, na.rm = T))
)

.filled.contour(
  x = lndscp_grid$x,
  y = lndscp_grid$y,
  z = lndscp_grid$z,
  levels = lndscp_levels,
  col    = topo.colors(length(lndscp_levels)-1)
)

coords <- fit$coords[chull(fit$coords),]
coords <- rbind(
  c(xlim[1], ylim[1]),
  c(xlim[1], ylim[2]),
  c(xlim[2], ylim[2]),
  c(xlim[2], ylim[1]),
  coords
)

dcoords <- geometry::delaunayn(coords)
dcoords <- dcoords[apply(dcoords, 1, function(x){ sum(x %in% 1:4) > 0 }),]

apply(dcoords, 1, function(x){
  polygon(
    x      = coords[,1][x],
    y      = coords[,2][x],
    col    = "grey90",
    border = "grey90"
  )
})

# polygon(fit$coords[chull(fit$coords),], border = "grey70", col = NA)

points(
  fit$coords,
  pch = 16,
  col = "grey80"
)
points(
  fit$coords
)

rect(
  xleft   = xlim[1],
  xright  = xlim[2],
  ybottom = ylim[1],
  ytop    = ylim[2],
  border  = "#000000",
  col     = NA,
  xpd     = TRUE
)
dev.off()

