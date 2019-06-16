
library(testthat)

## Test calculation of titer limits
testthat::context("Calculating titer limits")
pars <- ablandscape.control()

testthat::test_that("Limits with vectors", {
  
  # Check functioning of get_titer_lims with vectors
  expect_equal(object   = calc_titer_lims(titers = c("<10", 20, ">640")),
               expected = list(max_titers = c(-0.5, 1.5, pars$max.titer.possible),
                               min_titers = c(pars$min.titer.possible, 0.5, 6.5)))

})
    
testthat::test_that("Limits with matrices", {
  
  # Check functioning of get_titer_lims with matrices
  expect_equal(object   = calc_titer_lims(titers = matrix(c(10, 20, 40, 80, 160, 320), nrow = 2, ncol = 3)),
               expected = list(max_titers = matrix(0:5+0.5, nrow = 2, ncol = 3),
                               min_titers = matrix(0:5-0.5, nrow = 2, ncol = 3)))
  
})

testthat::test_that("Errors with a dataframe", {
  
  # Check errors of get_titer_lims
  expect_error(calc_titer_lims(titers = as.data.frame(matrix(c(10, 20, 40, 80, 160, 320), nrow = 2, ncol = 3))))
  
})
  
testthat::test_that("Limits with custom settings", {
  
  # Check setting of different defaults
  expect_equal(object   = calc_titer_lims(titers = c("<10", 20, ">640"), 
                                          control = list(min.titer.possible = -8,
                                                         max.titer.possible = 20)),
               expected = list(max_titers = c(-0.5, 1.5, 20),
                               min_titers = c(-8, 0.5, 6.5)))

})




## Test calculation of titer differences
testthat::context("Calculating titer differences")

testthat::test_that("Differences with vectors", {
  
  # Check functioning of get_titer_lims with vectors
  expect_equal(object   = calc_titer_diffs(titers1 = c("<10", 20, ">640"),
                                           titers2 = c(40, "<10", ">640")),
               expected = list(max_diffs = c(2.5 - pars$min.titer.possible, 
                                             -0.5 - 0.5, 
                                             pars$max.titer.possible - 6.5),
                               min_diffs = c(1.5 - -0.5, 
                                             pars$min.titer.possible - 1.5, 
                                             6.5 - pars$max.titer.possible)))
  
})

testthat::test_that("Differences with matrices", {
  
  # Check functioning of get_titer_lims with matrices
  expect_equal(object   = calc_titer_diffs(titers1 = matrix(c(10, 20, 40, 80, 160, 320), nrow = 2, ncol = 3),
                                           titers2 = matrix(c(10, 20, 40, 80, 160, 320), nrow = 2, ncol = 3)),
               expected = list(max_diffs = matrix(1, nrow = 2, ncol = 3),
                               min_diffs = matrix(-1, nrow = 2, ncol = 3)))
  
})

