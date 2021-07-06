test_that("Classify returns the expected factor", {
  # test on arbitrary thresholds:
  x <- simulated_umis %>%
    rule("T", "CD3E", ">", .1e-3) %>%
    rule("Treg", "FOXP3", ">", .01e-3, parent="T") %>%
    rule("Treg_act", "ICOS", ">", .01e-3, parent="Treg") %>%
    rule("B", "MS4A1", ">", .1e-3) %>%
    rule("B", "CD3E", "<", .1e-3)
  expect_equal(
    as.numeric(table(classify(obj=x, classes=c("T", "Treg", "Treg_act", "B")))),
    c(551, 0, 0, 649, 700))

})


test_that("classify's boolean output is as expected.", {
  obj <- simulated_umis %>% rule("T", "CD3E", ">", 1e-3) %>%
    rule("notT", "CD3E", "<", 1e-3) %>%
    rule("B", "MS4A1", ">", 1e-3, parent="notT") 
  res <- classify(obj, classes="B", return_logical_matrix=T )
  expect_equal(dim(res),
               dim(matrix(rep(FALSE,nrow(simulated_umis$raw)),
                          ncol=1)))
  expect_equal(sum(res), 145)
  expect_true(is.logical(res))
})




