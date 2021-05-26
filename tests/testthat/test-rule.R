obj_min <- list(raw=matrix(rpois(12, .5),ncol=3), embed=data.frame(u1=1:4, u2=4:1))


test_that("rule inputs are NULL in a healthy way", {
  null_error <- "rule does not accept NULL input for arguments: obj, class, feature and threshold."
  expect_error( rule(obj=NULL, class = "a", feature="b", threshold=1),      null_error )
  expect_error( rule(obj=obj_min, class = NULL, feature="b", threshold=1),  null_error )
  expect_error( rule(obj=obj_min, class = "a", feature=NULL, threshold=1),  null_error )
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=NULL),null_error )
  
})

test_that("rule inputs are of expected class", {
  # correct classes give no error:
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=1),  NA)
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



test_that("rule meets sanity checks", {
  expect_error( rule(obj=obj_min, class = "a", feature="b", t=.42, parent="a"),
                "Class and parent cannot be the same.")
  
})
