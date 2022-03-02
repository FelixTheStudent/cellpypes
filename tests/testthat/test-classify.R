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
    c(551, 26, 23, 649, 651))

})


test_that("replace_overlap_with may be one of the classes", {
  replace_with_class <- function(){
    simulated_umis %>%
      rule("CD3E+", "CD3E", ">", .5e-4) %>%
      rule("CD3E++","CD3E", ">", 5e-4) %>%
      classify(c("CD3E+","CD3E++"), replace_overlap_with="CD3E++")  
  }
  expect_error(replace_with_class(), NA)
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





test_that("overlap between class and its ancestry is not replaced.", {
  obj2 <- simulated_umis %>% 
    rule("B", "MS4A1", ">", 1e-4) %>% 
    rule("T", "CD3E", ">", .1e-4) %>% 
    rule("Treg", "FOXP3", ">", .05e-4, parent="T") 
    
  class_labels <- classify(obj2, classes=c("B","T","Treg"))
  # All Tregs overlap with Ts. In this case, I want Treg, not Unassigned:
  expect_true( any(class_labels =="Treg")) 
  # I solved bug triggered by 'class_res[,descendants, drop=F]', here's two tests for it:
  expect_error(classify(obj2, classes="T"), regexp = NA)  # test for error-free 
  expect_true(all(classify(obj2, classes="T") %in% c("T","Unassigned"))) 
  
  
  # When you implement common_parent, this error is supposed to remind you
  # to create same test as above, but with these class_labels:
  #     class_labels2 <- classify(obj2,
  #                               classes=c("B","T","Treg"), 
  #                               replace_overlap_with = "common_parent")
  expect_error(  classify(obj2, classes=c("B","T","Treg"), replace_overlap_with = "common_parent"),
                 regexp = "Not implemented yet, sorry.")
})


test_that("Non-existing parent gives intelligible error message.", {
  obj <- simulated_umis %>%
    rule("B", "MS4A1", ">", 1e-4) %>% 
    rule("T", "CD3E", ">", .1e-4) %>% 
    rule("Treg", "FOXP3", ">", .05e-4, parent="non-existing!") 
  expect_error(classify(obj), regexp = "A class has parent that does not exist -- double-check your rules!")
})
