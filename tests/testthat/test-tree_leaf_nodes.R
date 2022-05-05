test_that("tree_leaf_nodes returns expected leaves", {
  x <- simulated_umis
  x$raw <- rbind(x$raw, CD4=0, CXCL13=0, CD8B=0)
  x <- rule(x, "B", "MS4A1", ">", .1e-3) 
  x <- rule(x, "T", "CD3E", ">", .1e-3) 
  x <- rule(x, "Ttox", "CD8B",  ">", .1e-3, parent="T") 
  x <- rule(x, "T_CD4","CD4",   ">", .1e-3, parent="T") 
  x <- rule(x, "Treg", "FOXP3", ">", .1e-3, parent="T_CD4") 
  x <- rule(x, "Tfh",  "CXCL13",">", .1e-3, parent="T_CD4")
  expect_equal(
   tree_leaf_nodes(x$classes),
   c("B","Ttox", "Treg", "Tfh")
  )
})
