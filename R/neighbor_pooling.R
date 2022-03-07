

#' Title
#'
#' @param x 
#' @param neighbors 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom rlang is_double is_integer
pool_across_neighbors <- function(x, neighbors) {
  stopifnot(is_double(x) || is_integer(x))
  stopifnot(length(x)==nrow(neighbors))  
  
  if (ncol(neighbors)==nrow(neighbors)) {
    # assume neighbors is a graph matrix (typically with binary edge weights)
    return( as.numeric( as(neighbors, "dgCMatrix") %*% x ) )
  } else if (ncol(neighbors)<nrow(neighbors)) {
    # neighbors is assumed to be a matrix with neighbor indices for each cell
    rowSums(matrix(data=x[c(neighbors)], ncol=ncol(neighbors)))
  } else {
    stop("neighbors can't have more columns than rows.")
  }
}



#' Title
#'
#' @param obj 
#' @param class 
#' @param feature 
#' @param operator 
#' @param threshold 
#'
#' @return
#' @keywords internal
#'
#' @examples
evaluate_rule <- function(obj,
                          class,
                          feature,
                          operator,
                          threshold) {
  # check that obj has everything this rule needs:
  stopifnot(feature %in% rownames(obj$raw))
  
  # check_obj(obj) is called in functions calling evaluate_rule.

  K <- pool_across_neighbors(obj$raw[feature,], 
                             obj$neighbors)
  if (is.null(obj$totalUMI)) { 
    obj$totalUMI <- Matrix::colSums(obj$raw)
  } 
  S <- pool_across_neighbors(obj$totalUMI,
                             obj$neighbors)
    

  cdf <- ppois(K, S*threshold)
  switch(operator,
         ">" = cdf > .99,
         ">=" =cdf > .01,
         "<"  =cdf < .01,
         "<=" =cdf < .99)

}