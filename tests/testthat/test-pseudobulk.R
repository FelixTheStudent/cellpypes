
# Dummy celltypes and patients for simulated_umi, for testing:
df_meta <- data.frame(
  celltype = rep(c("X+Y-", "X+Y+", "X-Y+", "X-Y-"),
                 each = nrow(simulated_umis$embed)/4), # 4 cell types
  patient  = c("3", "500.", "*5", "/")
)

test_that("class_to_deseq2 example working", {
  expect_error( # NA means no error
    {
      data("simulated_umis") 
      # Meta data
      ncells <- ncol(simulated_umis$raw)
      dummy_variable <- function(x) factor(sample(x, ncells, replace=TRUE))
      meta_data <- data.frame(patient=dummy_variable(paste0("patient", 1:6)),
                              treatment=dummy_variable(c("control", "treated")))
      
      obj <- rule(simulated_umis, "T", "CD3E",">", 1)
      dds <- class_to_deseq2(obj, meta_data, "T", ~ treatment)
    }, NA)
})






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









# simulate dummy meta data for simulated_umis:
ncells <- ncol(simulated_umis$raw)
dummy_variable <- function(x) factor(sample(x, ncells, replace=TRUE))
meta_data <- data.frame(patient=dummy_variable(paste0("patient", 1:6)),
                        treatment=dummy_variable(c("control", "treated")))

test_that("class_to_deseq2 stops gracefully if class has zero cells.", {
  dds <- rule(simulated_umis, "T", "CD3E",">", 1000)
  expect_error(
    class_to_deseq2(dds, meta_data, "T", ~ treatment),
    "contains no cells"
    
  )
  
})


