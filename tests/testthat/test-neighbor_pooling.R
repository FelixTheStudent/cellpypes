set.seed(42)
dummy_dat <- matrix(rpois(3000, .1), ncol=30) # 30 cells without biological signal 
neighbors <- find_knn(dummy_dat, k = 15) # find 15 nearest neighbors 
neighbor_sum_gene1 <- pool_across_neighbors(dummy_dat[1,], neighbors$idx)

test_that("pool Docu example works", {
  expect_equal(
    neighbor_sum_gene1,
    c(3, 2, 1, 2, 1, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 2, 3, 2, 1, 2, 1, 2, 1, 1, 0,
      1, 0, 2, 2, 1)
  )
  
})




test_that("pooling with neighbor indices gives expected result", {
  expect_equal(41781, sum(pool_across_neighbors(simulated_umis$raw["CD3E",],
                                                simulated_umis$neighbors)))
  expect_equal( 2382, sum(pool_across_neighbors(simulated_umis$raw["FOXP3",],
                                                simulated_umis$neighbors)))
})






# simple nearest neighbor graph
square_matrix <- matrix(0, nrow = length(simulated_umis$celltype),
                        ncol = length(simulated_umis$celltype))
diag(square_matrix) <- 1
square_matrix[,1:5] <- 1

test_that("pooling with square matrix gives expected result", {
  expect_equal(16083, sum(pool_across_neighbors(simulated_umis$raw["CD3E",],
                                                square_matrix)))
  expect_equal( 104, sum(pool_across_neighbors(simulated_umis$raw["FOXP3",],
                                                square_matrix)))
})





test_that("evaluate_rule gives expected result", {
  expect_equal(602, sum(evaluate_rule(simulated_umis, "CD3E", ">", 1e-04)))
  expect_equal(64, sum(evaluate_rule(simulated_umis, "FOXP3", ">", 1e-04)))
})

