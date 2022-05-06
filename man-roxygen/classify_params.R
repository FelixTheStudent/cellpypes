#' @param classes Character vector with one or more class names.
#' If NULL (the default), plots finest available cell types
#' (all classes that are not parent of any other class).
#' @param replace_overlap_with Character string, by default: \code{"Unassigned"}.
#' See section \strong{Handling overlap}.
#' @param return_logical_matrix logical. If TRUE,
#' a logical matrix with
#' classes in columns and cells in rows is returned instead of resolving
#' overlaps with \code{replace_overlap_with}. 
#' If a single class is supplied, the matrix has exactly one
# column and the user can pipe it into "drop" to convert it to a vector.
