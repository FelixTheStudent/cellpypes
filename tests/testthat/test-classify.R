test_that("Classify returns the expected factor", {
  # test on arbitrary thresholds:
  x <- simulated_umis %>%
    rule("T", "CD3E", ">", .1e-3) %>%
    rule("Treg", "FOXP3", ">", .01e-3, parent="T") %>%
    rule("Treg_act", "ICOS", ">", .01e-3, parent="Treg") %>%
    rule("B", "MS4A1", ">", .1e-3) %>%
    rule("B", "CD3E", "<", .1e-3)
  expect_equal(
    classify(obj=x, classes=c("Treg_act", "B")),
    factor(c(rep("Treg_act",    23),
             rep("B",          649),
             rep("Unassigned",1228)
    ))
  )
})




