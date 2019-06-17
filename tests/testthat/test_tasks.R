
testthat::context("Tasks for today")

# Setup workspace
rm(list = ls())

# Load packages
library(Racmacs)
library(AbLandscapes)

# Read in a map
map <- h3basemap2015

# Fetch titer data
titers_pre  <- vaccine1997.titers.pre[12,]
titers_post <- vaccine1997.titers.post[12,]

# Fit a landscape
fit_pre <- ablandscape(titers    = titers_pre,
                       acmap     = map,
                       bandwidth = 10,
                       degree    = 1)

fit_post <- ablandscape(titers    = titers_post,
                        acmap     = map,
                        bandwidth = 10,
                        degree    = 1)

# View the fit
# view(fit)

# Get the summary path
summary_path <- read.csv("~/Desktop/snake_coords.csv", header = FALSE)
summary_path <- summary_path[seq(from = 1, to = nrow(summary_path), length.out = 100),]
summary_path

summary_lndscp_height_pre  <- predict(fit_pre,  coords = summary_path, crop2chull = FALSE)
summary_lndscp_height_post <- predict(fit_post, coords = summary_path, crop2chull = FALSE)

plot(x = seq_along(summary_lndscp_height_pre),
     y = summary_lndscp_height_pre)

plot(x = seq_along(summary_lndscp_height_post),
     y = summary_lndscp_height_post)

# view(fit)

