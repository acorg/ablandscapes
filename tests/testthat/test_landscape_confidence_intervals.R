
library(AbLandscapes)

## Test calculation of titer limits
testthat::context("Calculating confidence intervals")

warning("Need additional test for fitting confidence intervals")

# titers1 <- vaccine1997.titers.pre[10,]
# titers2 <- vaccine1997.titers.post[10,]
# coords  <- h3coords2015[match(names(titers2), rownames(h3coords2015)),]
# 
# # testthat::test_that("Predict new data with confidence intervals", {
# 
#   fit <- ablandscape.fit(titers     = titers2,
#                          coords     = coords,
#                          bandwidth  = 10,
#                          degree     = 1)
# 
#   output <- predict(fit, 
#                     coords = coords, 
#                     crop2chull = FALSE, 
#                     interval = "confidence",
#                     level = 0.95)
#   
#   testthat::expect_equal(ncol(output), 3)
# 
# # })


