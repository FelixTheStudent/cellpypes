#' @section knn_refine:
#' With `knn_refine > 0`, cellpypes refines cell type labels with a kNN classifier.
#' 
#' By default, cellpypes only assigns cells to a class if all relevant rules
#' apply.
#' In other words, all marker gene UMI counts in the cell's neighborhood all have to 
#' be clearly above/below their threshold.
#' Since UMI counts are sparse (even after neighbor pooling done by cellpypes),
#' this can leave many cells unassigned.
#' 
#' It is reasonable to assume an unassigned cell is of the same cell type as the
#' majority of its nearest neighbors.
#' Therefore, cellpypes implements a kNN classifier to further refine labels 
#' obtained by 
#' manually thresholding UMI counts.
#' `knn_refine = 0.3` means a cell is assigned the class label held by
#' most of its neighbors unless no class gets more than 30 %.
#' If most neighbors are unassigned, the cell will also be set to "Unassigned".
#' Choosing `knn_refine = 0.3` gives results reminiscent of clustering
#' (which assigns all cells),
#' while `knn_refine = 0.5` leaves cells 'in between' two similar
#' cell types unassigned.
#' 
#' 
#' We recommend looking at `knn_refine =  0` first as it's faster and
#' more directly tied to marker gene expression.
#' If assigning all cells is desired, we recommend `knn_refine = 0.3` or lower,
#' while `knn_refine = 0.5` makes cell types more 'crisp' by setting cells 
#' 'in between' related subtypes to "Unassigned".
