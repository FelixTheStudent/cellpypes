test_that("docu example works", {
  # Imagine we have 30 cells and 100 features:
  fmat <- matrix(rnorm(3000), ncol=30)
  nn <- find_knn(fmat,k=15)
  # nn$idx has 30 rows and 15 columns.
  expect_equal(dim(nn$idx), c(30, 15))
  expect_equal(dim(nn$dist), c(30, 15))
  
})
