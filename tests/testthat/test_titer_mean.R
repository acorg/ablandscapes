
## Test calculation of titer limits
testthat::context("Calculating maximum-likelihood means")

testthat::test_that("Mean titers", {
  
  logtiter.mean(c("40", "80", "<120"))
  
  testthat::expect_equal(object   = logtiter.mean(titers = c(20, 20, 20))$estimate$gmt,
                         expected = 1)
  
  testthat::expect_lt(object   = logtiter.mean(titers = "<10")$estimate$gmt,
                      expected = -1)
  
  testthat::expect_lt(object   = logtiter.mean(titers = c("<10", "<10", "<10"))$estimate$gmt,
                      expected = -1)
  
  testthat::expect_silent(object = logtiter.mean(titers = c("640", "1280", "1280"))$estimate$gmt)
  
})


testthat::test_that("Mean fold change", {
  
  logfoldchange.mean(titers1 = c("20", "80", "120"), 
                     titers2 = c("40", "80", "<120"))
  
  testthat::expect_equal(object   = logfoldchange.mean(titers1 = c(20, 20, 20),
                                                       titers2 = c(20, 20, 20))$estimate$gmt,
                         expected = 0)
  
})

