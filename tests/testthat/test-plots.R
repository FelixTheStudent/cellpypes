




test_that("plot_last does not require totalUMI", {
  obj_without_totalUMI <- list(
    raw=t(data.frame(CD3E=rpois(20, 5))),
    embed=data.frame(u1=1:20, u2=20:1),
    neighbors  =matrix(1:20, nrow=20, ncol=10))
  plot_WO_feat <- function() plot_last(
    obj_without_totalUMI %>% rule("T","CD3E",">",.1e-3),
    show_feat = FALSE)
  plot_W_feat <- function() plot_last(
    obj_without_totalUMI %>% rule("T","CD3E",">",.1e-3),
    show_feat = TRUE)
  expect_error(plot_WO_feat(), NA) # NA checks that there is no error
  expect_error(plot_W_feat(),  NA) # NA checks that there is no error
})



test_that("plot_classes handles embed tibbles", {
  obj <- simulated_umis
  obj$embed <- tibble::as_tibble(obj$embed)
  obj <- obj %>% rule("T", "CD3E", ">", 1e-4)
  expect_error(print(plot_classes(obj)), NA)
})

test_that(
  paste0("classify returns factor; otherwise ",
         "plot_classes has to use unique instead of levels"),
  {
    labels <- classify(rule(simulated_umis, "T", "CD3E", ">", .1e-3)) 
    expect_true(inherits(labels, "factor"))
  })



test_that("feat gives intelligible error messages", {
  res <- evaluate_promise(feat(simulated_umis, "CD3E", "MS4A1"))
  expect_true(grepl("Make sure to pass features as vector", res$warnings[1]))
  expect_true(grepl("Cannot convert object of class character into a grob", res$warnings[2]))
  
 
})
