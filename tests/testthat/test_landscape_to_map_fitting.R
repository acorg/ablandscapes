
library(AbLandscapes)

## Test calculation of titer limits
testthat::context("Fitting antibody landscapes to a map object")

titers1 <- vaccine1997.titers.pre[10,]
titers2 <- vaccine1997.titers.post[10,]
map     <- h3basemap2015

testthat::test_that("Fit landscape to map", {

  fit1 <- ablandscape(titers     = titers2,
                      acmap      = map,
                      bandwidth  = 10,
                      degree     = 1)
  
  testthat::expect_equal(
    sum(is.na(fit1$landscape$z)),
    0
  )
  
})
  
#   fit2 <- ablandscape.fit(titers     = titers2,
#                           coords     = coords,
#                           bandwidth  = 8,
#                           degree     = 1)
#   
#   fit3 <- ablandscape.fit(titers     = titers2,
#                           coords     = coords,
#                           bandwidth  = 10,
#                           degree     = 2)
# 
#   lessthans <- substr(titers2, 1, 1) == "<"
#   morethans <- substr(titers2, 1, 1) == ">"
# 
#   testthat::expect_equal(unique(fit1$residuals[lessthans | morethans]), NA_real_)
#   testthat::expect_equal(sum(is.na(fit1$residuals[!lessthans & !morethans])), 0)
#   testthat::expect_equal(sum(is.na(fit1$fitted.values)), 0)
#   testthat::expect_lt(sum(fit2$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
#   testthat::expect_lt(sum(fit3$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
# 
# })
# 
# 
# testthat::test_that("Fit delta landscape to coords", {
#   
#   fit1 <- ablandscape.delta.fit(titers1    = titers1,
#                                 titers2    = titers2,
#                                 coords     = coords,
#                                 bandwidth  = 10,
#                                 degree     = 1)
#   
#   fit2 <- ablandscape.delta.fit(titers1    = titers1,
#                                 titers2    = titers2,
#                                 coords     = coords,
#                                 bandwidth  = 8,
#                                 degree     = 1)
#   
#   fit3 <- ablandscape.delta.fit(titers1    = titers1,
#                                 titers2    = titers2,
#                                 coords     = coords,
#                                 bandwidth  = 10,
#                                 degree     = 2)
#   
#   measurable <- substr(titers1, 1, 1) != "<" &
#                 substr(titers2, 1, 1) != "<" &
#                 substr(titers1, 1, 1) != ">" &
#                 substr(titers2, 1, 1) != ">"
#   
#   testthat::expect_equal(unique(fit1$residuals[!measurable]), NA_real_)
#   testthat::expect_equal(sum(is.na(fit1$residuals[measurable])), 0)
#   testthat::expect_equal(sum(is.na(fit1$fitted.values)), 0)
#   testthat::expect_lt(sum(fit2$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
#   testthat::expect_lt(sum(fit3$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
#   
# })
# 
# 
# testthat::test_that("Predict new data", {
# 
#   fit <- ablandscape.fit(titers     = titers2,
#                          coords     = coords,
#                          bandwidth  = 10,
#                          degree     = 1)
# 
#   output <- predict(fit, coords = coords, crop2chull = FALSE)
#   testthat::expect_equal(fit$fitted.values[!is.na(fit$fitted.values)], output[!is.na(fit$fitted.values)])
# 
# })
# 
# 
# testthat::test_that("Predict new delta data", {
#   
#   fit <- ablandscape.delta.fit(titers1    = titers1,
#                                titers2    = titers2,
#                                coords     = coords,
#                                bandwidth  = 10,
#                                degree     = 1)
#   
#   output <- predict(fit, coords = coords, crop2chull = FALSE)
#   testthat::expect_equal(fit$fitted.values[!is.na(fit$fitted.values)], output[!is.na(fit$fitted.values)])
#   
# })

