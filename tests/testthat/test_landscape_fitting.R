
library(ablandscapes)
library(testthat)

## Test calculation of titer limits
context("Fitting antibody landscapes")

warning("Need test for when antigen coordinates has no rows")

titers1 <- vaccine1997.titers.pre[10,]
titers2 <- vaccine1997.titers.post[10,]
coords  <- h3coords2015[match(names(titers2), rownames(h3coords2015)),]

test_that("Fit landscape to coords", {

  fit1 <- ablandscape.fit(titers     = titers2,
                          coords     = coords,
                          bandwidth  = 10,
                          degree     = 1,
                          error.sd   = 1.1)
  
  fit2 <- ablandscape.fit(titers     = titers2,
                          coords     = coords,
                          bandwidth  = 8,
                          degree     = 1,
                          error.sd   = 1.1)
  
  fit3 <- ablandscape.fit(titers     = titers2,
                          coords     = coords,
                          bandwidth  = 10,
                          degree     = 2,
                          error.sd   = 1.1)

  lessthans <- substr(titers2, 1, 1) == "<"
  morethans <- substr(titers2, 1, 1) == ">"

  expect_equal(unique(fit1$residuals[lessthans | morethans]), NA_real_)
  expect_equal(sum(is.na(fit1$residuals[!lessthans & !morethans])), 0)
  expect_equal(sum(is.na(fit1$fitted.values)), 0)
  expect_lt(sum(fit2$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
  expect_lt(sum(fit3$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))

})


# test_that("Fit delta landscape to coords", {
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
#   expect_equal(unique(fit1$residuals[!measurable]), NA_real_)
#   expect_equal(sum(is.na(fit1$residuals[measurable])), 0)
#   expect_equal(sum(is.na(fit1$fitted.values)), 0)
#   expect_lt(sum(fit2$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
#   expect_lt(sum(fit3$residuals, na.rm = T), sum(fit1$residuals, na.rm = T))
#   
# })

warning("Need test for when prediction coordinates are given as a vector")

test_that("Predict new data", {

  fit <- ablandscape.fit(titers     = titers2,
                         coords     = coords,
                         bandwidth  = 10,
                         degree     = 1,
                         error.sd   = 1.1)

  output <- predict(fit, coords = coords, crop2chull = FALSE)
  expect_equal(fit$fitted.values[!is.na(fit$fitted.values)], output[!is.na(fit$fitted.values)])

})


test_that("Predict new delta data", {
  
  fit <- ablandscape.delta.fit(titers1    = titers1,
                               titers2    = titers2,
                               coords     = coords,
                               bandwidth  = 10,
                               degree     = 1,
                               error.sd   = 1.1)
  
  output <- predict(fit, coords = coords, crop2chull = FALSE)
  expect_equal(fit$fitted.values[!is.na(fit$fitted.values)], output[!is.na(fit$fitted.values)])
  
})


test_that("No warnings for >=", {
  
  expect_silent({
    ablandscape.fit(titers     = unlist(hanam.titers[hanam.info$subjectnumber == 17 & hanam.info$sampleyear == 2007,,drop=F]),
                    coords     = h3basemap2015$ag_coords[match(colnames(hanam.titers), h3basemap2015$ag_names),],
                    bandwidth  = 10,
                    degree     = 1,
                    error.sd   = 1.1)
  })
  
})


test_that("Predict grid data", {
  
  fit <- ablandscape.fit(
    titers     = titers2,
    coords     = coords,
    bandwidth  = 16,
    degree     = 2,
    error.sd   = 1.1
  )
  
  lndscp_grid <- predict_lndscp_grid(
    fit,
    crop2chull = FALSE
  )
  
  xlim <- range(lndscp_grid$x)
  ylim <- range(lndscp_grid$y)
  
})



