

#' Find approximate k-nearest neighbors
#'
#' @param featureMatrix Numeric matrix with features in rows, cells in columns.
#' Rows could be normalized genes or latent dimensions such as principal 
#' components.
#' @param k Number of neighbors to find.
#' @param n_trees RccpAnnoy builds a forest  of \code{n_trees} trees.
#' More trees gives higher precision when querying. Default: 50.
#' @param seed Random seed for neighbor search, default: 42.
#' 
#' 
#' @description Implements RcppAnnoy's approximate nearest neighbor search
#' (much faster than precise neighbors).
#' Random search is made reproducible using `set.seed(seed)`.
#' Hint: If you pass `find_knn`'s output directly to `uwot::umap` via the 
#' `nn_method` argument, make sure to set `umap`'s argument `n_sgd_threads`
#' to <=1 to ensure the UMAP embedding is reproducible.
#'
#' @return List with two slots: 
#' \itemize{
#'   \item \code{idx} A NxK matrix (N cells, K neighbors) containing the integer
#'   indexes of the approximate nearest neighbors in featureMatrix.
#'   Each cell is considered to be its own nearest neighbor, next to
#'   K-1 other neighbors.
#'   \item \code{dist} A NxK matrix containing the distances of the nearest neighbors.
#' }
#' Inspired by \code{uwot::umap}'s return value when setting \code{ret_nn=TRUE}.
#' 
#' @export
#'
#' @examples
#'   # Imagine we have 30 cells and 100 features:
#'  fmat <- matrix(rnorm(3000), ncol=30)
#'  nn <- find_knn(fmat,k=15)
#'  # nn$idx has 30 rows and 15 columns.
find_knn <- function(featureMatrix,
                     k=50,
                     n_trees = 50,
                     seed=42) {
  
  stopifnot(requireNamespace("RcppAnnoy", quietly = TRUE)) 
  
  # code below finds kNN between rows, not columns:
  featureMatrix <- t(featureMatrix)
  
  # Find nearest neighbors for UMAP and Louvain clustering:
  set.seed(seed) # seed ensures that UMAP gives reproducible result
  k_nn <- k 
  annoy <- methods::new( RcppAnnoy::AnnoyEuclidean, ncol(featureMatrix) )
  for( i in 1:nrow(featureMatrix) )
    annoy$addItem( i-1, featureMatrix[i,] )
  annoy$build( n_trees ) # builds a forest  of n_trees trees. More trees gives higher precision when querying.
  nn_cells <- t( sapply( 1:annoy$getNItems(), function(i) annoy$getNNsByItem( i-1, k_nn) + 1 ) )
  nndists_cells <- sapply( 1:ncol(nn_cells), function(j) sqrt( rowSums( ( featureMatrix - featureMatrix[ nn_cells[,j], ] )^2 ) ) )
  rm(annoy)
  
  return(list(idx=nn_cells,
              dist=nndists_cells))
}