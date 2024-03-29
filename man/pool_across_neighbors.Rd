% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighbor_pooling.R
\name{pool_across_neighbors}
\alias{pool_across_neighbors}
\title{Sum up x across neighbors in a nearest neighbor graph.}
\usage{
pool_across_neighbors(x, neighbors)
}
\arguments{
\item{x}{Numeric vector.}

\item{neighbors}{Nearest neighbor graph provided as NxK index matrix
(N observations, K neighbors) or NxN adjacency matrix.
Index matrices can be obtained with \link[cellpypes]{find_knn}
(specifically the slot \code{idx} in the list it returns).}
}
\value{
Numeric vector of length x.
}
\description{
Neighbor pooling means that x is summed across
the nearest neighbors.
}
\examples{
set.seed(42)
# simulate 30 cells without biological signal:
dummy_dat <- matrix(rpois(3000, .1), ncol=30) 
# find 15 approximate nearest neighbors 
neighbors <- find_knn(dummy_dat, k = 15) 
# pool gene1 counts across neighbors:
neighbor_sum_gene1 <- pool_across_neighbors(dummy_dat[1,], neighbors$idx)

}
