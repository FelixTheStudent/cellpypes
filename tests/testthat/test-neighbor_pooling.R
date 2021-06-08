
test_that("pooling with neighbor indices gives expected result", {
  expect_equal(41372, sum(pool_across_neighbors(simulated_umis$raw[,"CD3E"],
                                                simulated_umis$neighbors)))
  expect_equal( 2193, sum(pool_across_neighbors(simulated_umis$raw[,"FOXP3"],
                                                simulated_umis$neighbors)))
})






# simple nearest neighbor graph
square_matrix <- matrix(0, nrow = length(simulated_umis$celltype),
                        ncol = length(simulated_umis$celltype))
diag(square_matrix) <- 1
square_matrix[,1:5] <- 1

test_that("pooling with square matrix gives expected result", {
  expect_equal(16083, sum(pool_across_neighbors(simulated_umis$raw[,"CD3E"],
                                                square_matrix)))
  expect_equal( 104, sum(pool_across_neighbors(simulated_umis$raw[,"FOXP3"],
                                                square_matrix)))
})
