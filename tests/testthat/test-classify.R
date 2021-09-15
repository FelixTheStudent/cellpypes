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




test_that("classify does not require totalUMI", {
  obj_without_totalUMI <- list(
    raw=data.frame(CD3E=rep(2, 20)),
    embed=data.frame(u1=1:20, u2=20:1),
    neighbors  =matrix(1:20, nrow=20, ncol=10))
  # NA checks that there is no error:
  expect_error(classify(obj_without_totalUMI %>% rule("T","CD3E",">",.1e-3)),
               NA)
  # Next to running without error, classify should return class labels:
  expect_equal( as.character(classify(obj_without_totalUMI %>% rule("T","CD3E",">",.1e-3))),
                rep("T", 20))
})
