# Functions for forming pseudobulks ("ensembls") by summing up single cells





# Internal helper function. 
# 
# I had this function for years, so I wanted to keep it instead of incorporating
# it into the \code{pseudobulk} function.
#
# @param factor A factor you want to use for pseudobulking.
#
# @return A matrix where each factor level has been turned into a binary
# column. 
fac2matrix <- function(factor){
    tmp <- factor
    p <- stats::model.matrix(~0+tmp)
    colnames(p) <- gsub("^tmp", "", colnames(p))
    return(p)
}
 
  
  
#' Form pseudobulks from single cells.
#'
#' @param raw A matrix with raw UMI counts, cells in columns. 
#' @param pseudobulk_id A factor that identifies which cells should go to
#' which pseudobulk. Generate pseudobulk_ids with the \link[cellpypes]{pseudobulk_id} function! 
#' 
#' @description 
#' Sum up cells in count matrix \code{raw} for bulk RNA methods such as
#' DESeq2.
#'
#' @return A matrix where each column is a pseudobulk and each row a gene.
#' @export
#'
#' @examples
#' # Create pseudobulk counts and coldata for DESeq2: 
#' coldata <- data.frame(
#'   celltype = rep(c("X+Y-", "X+Y+", "X-Y+", "X-Y-"),
#'                  each = nrow(simulated_umis$embed)/4), # 4 cell types
#'   patient  = c("3", "500.", "*5", "/")
#' )
#' coldata$pseudobulk_id <- pseudobulk_id(coldata)
#' counts <- pseudobulk(simulated_umis$raw, coldata$pseudobulk_id)
#' # Use counts/coldata as input for DESeqDataSetFromMatrix (DESeq2).
#' 
#' @importMethodsFrom Matrix %*%
pseudobulk <- function(raw, pseudobulk_id) {
  raw %*% fac2matrix( pseudobulk_id )
}








#' Generate unique IDs to identify your pseudobulks.
#'
#' @param factor_df Data frame where each column helps to identify a pseudobulk.
#' Each row in factor_df corresponds to a single cell in your raw count matrix.
#' Typical columns of factor_df are for example patient, treatment and cell type -- anything
#' that uniquely identifies a replicate.
#' @return Factor with syntactically valid and unique IDs.
#' @details Wraps \link[base]{make.names} to generate syntactically valid IDs.
#' Use these IDs in the \link[cellpypes]{pseudobulk} function.
#' Note that this function combines all columns in factor_df, so only include
#' the columns that uniquely identify replicates.
#' Cells from the same experimental unit
#  (combination of patient, treatment, celltype, etc.) 
#  get the same pseudobulk_id.
#' 
#' @description 
#' This function generates unique IDs that are valid colnames as well.
#' Use these IDs in function pseudobulk.
#' 
#' @export
#'
#' @examples
#' # Create pseudobulk counts and coldata for DESeq2: 
#' coldata <- data.frame(
#'   celltype = rep(c("X+Y-", "X+Y+", "X-Y+", "X-Y-"),
#'                  each = nrow(simulated_umis$embed)/4), # 4 cell types
#'   patient  = c("3", "500.", "*5", "/")
#' )
#' coldata$pseudobulk_id <- pseudobulk_id(coldata)
#' counts <- pseudobulk(simulated_umis$raw, coldata$pseudobulk_id)
#' # Use counts/coldata as input for DESeqDataSetFromMatrix (DESeq2).
pseudobulk_id <- function(factor_df) {
 
  id <- do.call(paste, factor_df)
  # machine-readable and unique:
  id_df <- data.frame(id = unique(id),
                      pseudobulk_id = make.names(unique(id), unique=TRUE))
  # pseudobulk_ids fit factor_df (which may have multiple cells with same ID):
  colnames(factor_df)[colnames(factor_df)=="id"] <- "dummy_name2"
  factor_df <- cbind(factor_df, id=id)
  # preserving order after merge requires this:
  colnames(factor_df)[colnames(factor_df)=="rank"] <- "dummy_name"
  factor_df$rank <- 1:nrow(factor_df)
  res <- merge(factor_df, id_df, by="id", sort=FALSE)
  # merge does not presever order, let's make sure:
  res <- res[order(res$rank),]
  
  factor(res$pseudobulk_id)
  
}






#' Create DESeq2 object for a given cell type
#'
#' @template param_obj
#' @param meta_df Data frame where each column helps to identify a pseudobulk.
#' Typical columns of meta_df are for example patient, treatment and
#' cell type -- anything
#' that uniquely identifies a replicate / batch / 10x run. 
#' Each row in meta_df corresponds to a single cell in your raw count matrix.
#' @param class The name of cellpypes class for which you want to test 
#' for differential expression.
#' @param design A formula based on columns in \code{meta_df}.
#' To test differential expression between two groups
#' in meta_df$condition, use formula \code{~ condition}.
#' More complex formulas (e.g. with interactions) are possible, for example 
#' \code{~ genotype + treatment + genotype:treatment}.
#' 
#' @description 
#' Create a DESeq2 data set (`dds' in the 
#' \href{https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html}{DESeq2 vignette})
#' for the specified class (cell type).
#'
#' @return A DESeq2 object (e.g. dds)
#' @template cellpypes_obj
#' @export
#' 
#' @importFrom rlang is_scalar_character
#'
#' @examples
#' data("simulated_umis") 
#' # Meta data
#' ncells <- ncol(simulated_umis$raw)
#' dummy_variable <- function(x) factor(sample(x, ncells, replace=TRUE))
#' meta_data <- data.frame(patient=dummy_variable(paste0("patient", 1:6)),
#'                         treatment=dummy_variable(c("control", "treated")))
#' 
#' obj <- rule(simulated_umis, "T", "CD3E",">", 1e-4)
#' \donttest{ # > 5 s in CRAN check
#' dds <- class_to_deseq2(obj, meta_data, "T", ~ treatment)
#' }
class_to_deseq2 <- function(obj, meta_df, class, design = ~ condition) {
  deseq_status <- requireNamespace("DESeq2", quietly = TRUE)
  if(!deseq_status) stop("Install DESeq2 to use class_to_deseq2 function.")
  check_obj(obj)
  stopifnot(
    inherits(meta_df, "data.frame"),
    inherits(design, "formula"),
    is_scalar_character(class),
    any(obj$classes$class == class))
  
  # After the fact, colnames(dds) should reveal which cell type was used:
  meta_df$celltype <- class
  
  # dds should only contain the relevant cell type to keep it simple:
  is_class <- classify(obj, classes=class) == class
  if(sum(is_class)==0) stop(paste0("Cell type '",class,"' contains no cells."))
  bulks <- pseudobulk(obj$raw[,is_class], pseudobulk_id(meta_df[is_class,]))
  # I don't want to depend on dplyr (huge package), so I write 'distinct' myself:
  coldat <- meta_df[is_class,]
  coldat$rowID <- pseudobulk_id(meta_df[is_class,])
  coldat <- coldat[!duplicated(coldat$rowID),] 
  rownames(coldat) <- coldat$rowID
  coldat <- coldat[,colnames(coldat) != "rowID", drop=FALSE]
  # bring to same order (good practice and required for DESeq2)
  coldat <- coldat[colnames(bulks), , drop=FALSE] 
  
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = bulks,
    colData = coldat,
    design = design)
  dds
}







