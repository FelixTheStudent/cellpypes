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
#' @template cellpypes_obj
#' @export
pype_from_seurat <- function(seurat) {
  seurat_status <- requireNamespace("Seurat", quietly=TRUE) &&
    requireNamespace("SeuratObject", quietly = TRUE)
  if(!seurat_status) stop("Install Seurat to use this function.")
  stopifnot(inherits(seurat, "Seurat"))
  
  
  # I pick graphs according to this wish list (order matters):
  # WNN_snn, SCT_snn, RNA_snn.
  all_graphs <- names(seurat@graphs)
  if (is.null(all_graphs)) {
    seurat <- Seurat::FindNeighbors(seurat)
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
  } else if (any(grepl("umap|tsne", names(seurat@reductions)))) {
    # reductions may have arbitrary names, e.g. rna50.umap
    use <- names(seurat@reductions)[grepl("umap|tsne|UMAP|TSNE", names(seurat@reductions))]
    use <- use[1] # in case there are multiple reductions
    dimension_names <- colnames(seurat@reductions[[use]]@cell.embeddings)
  } else {
    stop("Neither UMAP nor tSNE found.")
  }
  
  list(
    raw      =SeuratObject::GetAssayData(seurat, "counts"),
    neighbors=methods::as(seurat@graphs[[graph_choice]], "dgCMatrix")>.1,
    embed    =Seurat::FetchData(seurat, dimension_names),
    totalUMI = seurat$nCount_RNA
  ) 
} 