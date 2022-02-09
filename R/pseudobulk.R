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
    p <- model.matrix(~0+tmp)
    colnames(p) <- gsub("^tmp", "", colnames(p))
    return(p)
}
 
  
  
#' Form pseudobulks from single cells.
#'
#' @param raw A matrix with raw UMI counts, cells in columns. 
#' @param pseudobulk_id A factor that identifies which cells should go to
#' which pseudobulk. Generate pseudobulk_ids with the \link[cellpypes]{pseudobulk_id} function! 
#'
#' @return
#' @export
#'
#' @examples
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








#' Generate syntactically valid and unique IDs to identify your pseudobulks.
#'
#' @param factor_df Data frame where each column helps to identify a pseudobulk.
#' Each row in factor_df usually corresponds to a single cell in your raw count matrix.
#' Typical columns of factor_df are for example patient, treatment and cell type -- anything
#' that uniquely identifies a replicate.
#' @return Factor with syntactically valid and unique IDs.
#' @details Wraps \link[base]{make.names} to generate syntactically valid IDs.
#' Use these IDs in the \link[cellpypes]{pseudobulk} function.
#' Note that this function combines all columns in factor_df, so only include
#' The columns that uniquely identify replicates.
#' Cells from the same experimental unit
#  (combination of patient, treatment, celltype, etc.) 
#  get the same pseudobulk_id.
#' 
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










