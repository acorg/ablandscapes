
## Test calculation of titer limits
testthat::context("Calculating maximum-likelihood means")

testthat::test_that("Mean titers", {
  
  testthat::expect_equal(object   = mean.titer(titers = c(20, 20, 20))$estimate$gmt,
                         expected = 1)
  
})


testthat::test_that("Mean fold change", {
  
  testthat::expect_equal(object   = mean.foldchange(titers1 = c(20, 20, 20),
                                                    titers2 = c(20, 20, 20))$estimate$gmt,
                         expected = 0)
  
})

