# Idea behind this function:
# Converting Seurat objects to cellpypes object is best done with a wrapper 
# that guides the user towards reasonable choices.




#' Convert Seurat to cellpypes object. 
#' 
#' Start cellpyping a Seurat object.
#' This function saves the user from building his own cellpypes object,
#' which is done with \code{list(umi, neighbors,embed, totalUMI)}.
#'   
#'
#' @param seurat A Seurat object.
#'
#' @return A cellpypes object.
#' @export
#'
#' @examples
pype_from_seurat <- function(seurat) {
  stopifnot(inherits(malt, "Seurat"))
  
  
  # I pick graphs according to this wish list (order matters):
  # WNN_snn, SCT_snn, RNA_snn.
  all_graphs <- names(seurat@graphs)
  if (is.null(all_graphs)) {
    seurat <- FindNeighbors(seurat)
    all_graphs <- names(seurat@graphs)
  }
  snn_graphs <- all_graphs[grepl("_snn$", all_graphs)]
  graph_choices <- paste(c("WNN", "SCT", "RNA"), "snn", sep="_")
  graph_choice <- graph_choices[graph_choices %in% snn_graphs][1]
  
  # pick embedding according to this wish list order: umap, then tsne
  if ("umap" %in% names(seurat@reductions)) {
    dimension_names <- c("UMAP_1", "UMAP_2")
  } else if ("tsne" %in% names(seurat@reductions)) {
    dimension_names <- c("tSNE_1", "tSNE_2")
  } else {
    stop("Neither UMAP nor tSNE found.")
  }
  
  list(
    raw      =Matrix::t(GetAssayData(seurat, "counts")),
    neighbors=as(seurat@graphs[[graph_choice]], "dgCMatrix")>.1,
    embed    =FetchData(seurat, dimension_names),
    totalUMI = seurat$nCount_RNA
  ) 
} 