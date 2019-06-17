
rm(list = ls())

## Test calculation of titer limits
testthat::context("Plotting antibody landscapes")

titers1 <- vaccine1997.titers.pre[10,]
titers2 <- vaccine1997.titers.post[10,]
coords  <- h3coords2015[match(names(titers2), rownames(h3coords2015)),]

lndscp <- ablandscape(titers    = titers2,
                      acmap     = h3basemap2015,
                      bandwidth = 10,
                      degree    = 1)

data3js <- AbLandscapes::view(lndscp)

r3js::debug3js(data3js, filename = "ablandscape.js")
