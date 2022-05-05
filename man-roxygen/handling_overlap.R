#' @section Handling overlap:
#' Overlap denotes all cells
#' for which rules from multiple classes apply, and these cells will be
#' labeled as \code{Unassigned} by default.
#' If you are in fact interested in where the overlap is,
#' set \code{return_logical_matrix}=TRUE and inspect the result.
#' Note that 
#' it matters whether you call \code{classify("Tcell")} or
#' \code{classify(c("Tcell","Bcell")} – any existing overlap between T and B cells
#' is labelled as \code{Unassigned} in 
#' this second call, but not in the first.
#' 
#' Replacing overlap happens only between mutually 
#' exclusive labels (such as Tcell and Bcell), but
#' not within a lineage.
#' To make an example, overlap is NOT replaced between child (PD1+Ttox) and
#' parent (Ttox) or any other ancestor (Tcell), but instead the 
#' most detailed cell type (PD1+Ttox) is returned.
#' 
#' All of the above is also true for \code{plot_classes}, as it wraps \code{classify}.

