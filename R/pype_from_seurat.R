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
#' @param graph_name Supply one of the graphs. To see options, type 
#' `names(seurat@graphs)`. If left empty (`NULL`, the default),
#' `pype_from_seurat` will try to guess the correct name for you.
#'
#' @return A cellpypes object.
#' @template cellpypes_obj
#' @export
pype_from_seurat <- function(seurat, graph_name=NULL) {
  seurat_status <- requireNamespace("Seurat", quietly=TRUE) &&
    requireNamespace("SeuratObject", quietly = TRUE)
  if(!seurat_status) stop("Install Seurat to use this function.")
  stopifnot(inherits(seurat, "Seurat"))
  graph_name_status = rlang::is_scalar_character(graph_name) || is.null(graph_name)
  stopifnot("graph_name must be a character scalar."=graph_name_status)
  
  if(is.null(graph_name)) {
    # I pick graphs according to this wish list (order matters):
    # WNN_snn, SCT_snn, RNA_snn.
    all_graphs <- names(seurat@graphs)
    if (is.null(all_graphs)) {
      seurat <- Seurat::FindNeighbors(seurat)
      all_graphs <- names(seurat@graphs)
    }
    
    # order matters:
    seurat_graph_prefixes <- c("integrated","WNN", "SCT", "RNA") 
    graph_choices <- c(paste0(seurat_graph_prefixes, "_snn"),
                       paste0(seurat_graph_prefixes, "_nn"))
    
    graph_choices <- graph_choices[graph_choices %in% names(seurat@graphs)]
    if(length(graph_choices)>0) {
      graph_choice <- graph_choices[1]
    } else {
      stop("No neighbor graph found. Try passing names(seurat@graphs) to argument graph_name.")
    }
    
  } else { # user-supplied graph_name
    stopifnot("graph_name not present in names(seurat@graphs)"=graph_name %in% names(seurat@graphs))
    graph_choice <- graph_name
  }
 
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
