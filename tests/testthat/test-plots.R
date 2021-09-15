




test_that("plot_last does not require totalUMI", {
  obj_without_totalUMI <- list(
  raw=data.frame(CD3E=rpois(20, 5)),
  embed=data.frame(u1=1:20, u2=20:1),
  neighbors  =matrix(1:20, nrow=20, ncol=10))
  # NA checks that there is no error:
  expect_error(plot_last(obj_without_totalUMI %>% rule("T","CD3E",">",.1e-3)),
               NA)
})


