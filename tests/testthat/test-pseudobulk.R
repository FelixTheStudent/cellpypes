
# Dummy celltypes and patients for simulated_umi, for testing:
df_meta <- data.frame(
  celltype = rep(c("X+Y-", "X+Y+", "X-Y+", "X-Y-"),
                 each = nrow(simulated_umis$embed)/4), # 4 cell types
  patient  = c("3", "500.", "*5", "/")
)








test_that("pseudobulk_ids handle special characters well", {
  # pseudobulk_ids should have 16 levels (4 celltype x 4 patient):
  expect_equal(length(levels( pseudobulk_id(df_meta) )), 16)
})


test_that("Example in pseudobulk_ids docu is working", {
  coldata <- df_meta
  expect_error( { # expect this runs without error
   coldata$pseudobulk_id <- pseudobulk_id(coldata)
   counts <- pseudobulk(simulated_umis$raw, coldata$pseudobulk_id)   
  }, NA)


})
