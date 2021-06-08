

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


# evaluate_rule <- function(obj,
#                           class,
#                           feature,
#                           operator,
#                           threshold) {
#  
# 
# rule <- sobj$rules[1,,drop=TRUE]
# K <- pool_nn(pcpc$raw[, rule$gene], pcpc$neighbors)
# S <- pool_nn(pcpc$totals,           pcpc$neighbors)
# cdf <- ppois(K, S*rule$threshold)
# switch(rule$operator,
#        ">" = cdf > .99,
#        ">=" =cdf > .01,
#        "<"  =cdf < .01,
#        "<=" =cdf < .99)
#   
# }