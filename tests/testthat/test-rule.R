
# example from docu (man file):
# T cells are CD3E+
obj <- rule(simulated_umis, "T", "CD3E", ">", .1)
# T cells are MS4A1-
obj <- rule(obj, "T", "MS4A1", "<", 1)
# Tregs are a subset of T cells:
obj <- rule(obj, "Treg", "FOXP3", ">", .1, parent="T") 

test_that("rule Docu example works", {
  expect_equal(tabulate(classify(obj, c("T", "Treg"))),
               c(445, 42, 1413))
})

test_that("mixing CP10K and fractions works", {
  obj2 <- rule(obj, "T", "MS4A1", "<", 1e-4, use_CP10K = FALSE)
  expect_equal(classify(obj),
               classify(obj2))
})



obj_min <- list(raw=t(matrix(rpois(12, .5),ncol=3)), embed=data.frame(u1=1:4, u2=4:1))


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
  class_error <- "^Wrong format"
  expect_error( rule(obj=obj_min, class = 42, feature="b", threshold=1), 
                "^class must be a single string")
  expect_error( rule(obj=obj_min, class = "a", feature=42, threshold=1),  
                "^feature must be a single string")
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=1, operator=42),  
                "^operator must be a single string")
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold="forty-two"),
                "^threshold must be a single number")
  expect_error( rule(obj=obj_min, class = "a", feature="b", threshold=42, parent=5),
                "^parent must be NULL or a single string")
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

test_that("rule can assign same feature to different classes.", {
  obj <- simulated_umis 
  obj <- rule(obj, "T", "CD3E", ">", 42)
  obj <- rule(obj, "B", "MS4A1",">", 42)
  obj <- rule(obj, "B", "CD3E", "<", 42)
  expect_equal(obj$rules,
               data.frame(
                 class=c("T","B","B"), 
                 feature= c("CD3E", "MS4A1", "CD3E"),
                 operator=c(">",">","<"),
                 threshold=42e-4 # e-4 makes it resemble internal CP10K value
               ))
})


test_that("rule moves modified rule to position of most recent rule.", {
  x <- simulated_umis 
  x <- rule(x, "T", "CD3E", ">",    42) 
  x <- rule(x, "B", "MS4A1",">",    42)
  x <- rule(x, "Treg", "FOXP3",">", 42, parent="T") 
  x <- rule(x, "T", "CD3E", "<",    41)
  expect_equal(x$rules$class,
               c("B", "Treg", "T")
  )
})

test_that("Existing classes and existing features are handled correctly", {
  hasT <- rule(obj=simulated_umis, class="T", feature="CD3E", operator="<",
               threshold=9001)
  hasT <- rule(obj=hasT, class="T", feature="CD3E", operator=">",
               threshold=42)
  expect_equal(hasT$classes, data.frame(class="T", parent="..root.."))
  expect_equal(hasT$rules, data.frame(class="T",
                                      feature="CD3E",
                                      operator=">",
                                      # e-4 makes it resemble internal CP10K value:
                                      threshold=42e-4)) 
  hasT <- rule(obj=hasT, class="T", feature="MS4A1", operator="<",
               threshold=42) 
  expect_equal(hasT$classes, data.frame(class="T", parent="..root.."))
  expect_equal(hasT$rules, data.frame(
    class=c("T","T"),
    feature=c("CD3E", "MS4A1"),
    operator=c(">", "<"),
    threshold=42e-4) ) # e-4 makes it resemble internal CP10K value
})
  

test_that("Existing feature is handled well even if obj is saved into variable.", {
  x <- simulated_umis
  x <- rule(x, "T", "CD3E", threshold = .1e-3) 
  x <- rule(x, "Treg", "FOXP3", threshold = .1e-3)
  x <- rule(x, "T", "CD3E", threshold = .1e-3)
  expect_true(is_rules(x$rules))
  
})



