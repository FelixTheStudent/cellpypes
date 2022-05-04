#' @section Handling overlap:
#' If cell types overlap, classify
#' returns \code{Unassigned} for cells with
#' mutually exclusive labels (such as Tcell and Bcell).
#' For this reason it matters whether you call \code{classify("Tcell")} or
#' \code{classify(c("Tcell","Bcell")} – any existing overlap between T and B cells
#' is labelled as \code{Unassigned} in 
#' this second call, but not in the first.
#' If a cell gets multiple labels but from the same lineage (e.g. Tcell and CD8+ T), the more detailed class is returned (CD8+ T).
#' All of the above is also true for \code{plot_classes}, as it wraps \code{classify}.
