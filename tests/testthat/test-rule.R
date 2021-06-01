obj_min <- list(raw=matrix(rpois(12, .5),ncol=3), embed=data.frame(u1=1:4, u2=4:1))


test_that("rule checks whether its inputs are NULL in a healthy way", {
  null_error <- "rule does not accept NULL input for arguments: obj, class, feature and threshold."
  expect_error( rule(obj=NULL, class = "a", feature="b", threshold=1),      null_error )
  expect_error( rule(obj=obj_min, class = NULL, feature="b", threshold=1),  null_error )
  expect_error( rule(obj=obj_min, class = "a", feature=NULL, threshold=1),  null_error )
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=NULL),null_error )
  
})

test_that("rule verifies its inputs are of expected class", {
  # correct classes give no error:
  #    I'm not testing this currently, because the obj structure is still changing.
  
  # a wrong class in the input arguments gives error:
  class_error <- "^all"
  expect_error( rule(obj=obj_min, class = 42, feature="b", threshold=1),  class_error )
  expect_error( rule(obj=obj_min, class = "a", feature=42, threshold=1),  class_error )
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=1, operator=42),  class_error )
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold="forty-two"),class_error )
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=42, parent=5),class_error )
  # currently I don't check obj:
  # expect_error( rule(obj="lol", class = "a", feature="b", threshold=1),      class_error )
})



test_that("rule has sanity checks in place.", {
  expect_error( rule(obj=obj_min, class = "a", feature="b", t=.42, parent="a"),
                "Class and parent cannot be the same.")
  
})



test_that("rule adds rules as intended", {
  hasT <- rule(obj=simulated_umis, class="T", feature="CD3E", operator=">", threshold=42)
  # class can't be added twice:
  hasT <- rule(obj=hasT,           class="T", feature="CD3E", operator=">", threshold=42)
  expect_equal(hasT$classes,
               data.frame(class="T", parent="..root.."))
  
  hasB <- rule(obj=hasT,           class="B", feature="MS4A1", operator="<", threshold=42)
  hasM <- rule(obj=hasB,           class="M", feature="MS4A1", operator="<", threshold=42)
  # parent can be changed in retrospect:
  hasM <- rule(obj=hasM, class="M", feature="MS4A1", operator="<", threshold=42,
               parent="B")
  expect_equal(hasM$classes, data.frame(class=c("T", "B", "M"), 
                                        parent=c("..root..","..root..", "B")))
})

test_that("Existing classes and existing features are handled correctly", {
  hasT <- rule(obj=simulated_umis, class="T", feature="CD3E", operator="<",
               threshold=9001)
  hasT <- rule(obj=hasT, class="T", feature="CD3E", operator=">",
               threshold=42)
  expect_equal(hasT$classes, data.frame(class="T", parent="..root.."))
  expect_equal(hasT$rules, data.frame(class="T", feature="CD3E", operator=">",
                                      threshold=42))
  hasT <- rule(obj=hasT, class="T", feature="MS4A1", operator="<",
               threshold=42) 
  expect_equal(hasT$classes, data.frame(class="T", parent="..root.."))
  expect_equal(hasT$rules, data.frame(
    class=c("T","T"), feature=c("CD3E", "MS4A1"), operator=c(">", "<"),threshold=42) )
})
  



