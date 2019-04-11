context("test-frequentdirections.R")

test_that("all_zero_row_index() works", {
  eps <- 10^(-3)
  x <- matrix(c(1,1,0,10^(-6),2,2), nrow = 3, byrow = TRUE)
  testthat::expect_equal(all_zero_row_index(x, 10^(-3)), 2)
  testthat::expect_equal(all_zero_row_index(x, 10^(-8)), integer(0))
})
