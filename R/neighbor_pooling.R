

#' Sum up x across neighbors in a nearest neighbor graph.
#'
#' @param x Numeric vector.
#' @param neighbors Nearest neighbor graph provided as NxK index matrix 
#' (N observations, K neighbors) or NxN adjacency matrix.
#' Index matrices can be obtained with \link[cellpypes]{find_knn}
#' (specifically the slot \code{idx} in the list it returns). 
#' 
#' @description Neighbor pooling means that x is summed across 
#' the nearest neighbors.
#'
#' @return Numeric vector of length x.
#' @export
#'
#' @examples
#' set.seed(42)
#' # simulate 30 cells without biological signal:
#' dummy_dat <- matrix(rpois(3000, .1), ncol=30) 
#' # find 15 approximate nearest neighbors 
#' neighbors <- find_knn(dummy_dat, k = 15) 
#' # pool gene1 counts across neighbors:
#' neighbor_sum_gene1 <- pool_across_neighbors(dummy_dat[1,], neighbors$idx)
#' 
#' @importFrom rlang is_double is_integer
pool_across_neighbors <- function(x, neighbors) {
  stopifnot("x must be numeric"=is.numeric(x))
  stopifnot("length of x must match nrow of neighbors"=length(x)==nrow(neighbors))  
  
  if (ncol(neighbors)==nrow(neighbors)) {
    # assume neighbors is a graph matrix (typically with binary edge weights)
    return( as.numeric( methods::as(neighbors, "dgCMatrix") %*% x ) )
  } else if (ncol(neighbors)<nrow(neighbors)) {
    # neighbors is assumed to be a matrix with neighbor indices for each cell
    rowSums(matrix(data=x[c(neighbors)], ncol=ncol(neighbors)))
  } else {
    stop("neighbors can't have more columns than rows.")
  }
}



#' Evaluate rule to obtain positive / negative cells 
#'
#' @template param_obj
#' @param feature Character scalar naming the gene you'd like to threshold. 
#' @template param_operator
#' @param threshold Numeric scalar with the expression threshold separating positive
#' from negative cells. Evaluate_rule expects a plain fraction, not CP10K
#' (while rule has the use_CP10K argument).
#' @param overdispersion Defaults to 0.01, only change it if you know
#' what you are doing.
#' If set to 0, the NB simplifies to the Poisson distribution, and larger
#' values give more variance.
#' The 0.01 default value follows the recommendation by
#' Lause, Berens and Kobak (Genome Biology 2021) to use
#' `size=100` in \link[stats]{pnbinom} for typical data sets.
#' 
#' 
#' @description The rule is defined by feature, operator and threshold.
#'
#' @return logical vector stating for each cell whether the rule is true.
#' 
#' @template cellpypes_obj
#' 
#' @keywords internal
#'
evaluate_rule <- function(obj,
                          feature,
                          operator,
                          threshold,
                          overdispersion=0.01) {
  # This is a separate function for two reasons:
  #    * count pooling with NB is just one of many options I might implement
  #    * it's being used by classify AND by plot_last, so separate function.
  
  
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
    
  # before cellpypes version 0.1.14:
  # cdf <- stats::ppois(K, S*threshold)
  cdf <- stats::pnbinom(K, mu=S*threshold, size=1/overdispersion)
  switch(operator,
         # >= and <= are currently prevented by stopifnot in rule
         ">" = cdf > .99,
         ">=" =cdf > .01,
         "<"  =cdf < .01,
         "<=" =cdf < .99)

}