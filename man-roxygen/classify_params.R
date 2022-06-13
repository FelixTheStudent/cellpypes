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
#' column and the user can pipe it into "drop" to convert it to a vector.
#' @param overdispersion Defaults to 0.01, only change it if you know
#' what you are doing.
#' If set to 0, the NB simplifies to the Poisson distribution, and larger
#' values give more variance.
#' The 0.01 default value follows the recommendation by
#' Lause, Berens and Kobak (Genome Biology 2021) to use
#' `size=100` in \link[stats]{pnbinom} for typical data sets.
