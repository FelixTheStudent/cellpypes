#' @section cellpypes Objects:
#' A cellpypes object is a \link[base]{list} with four slots:
#' \describe{
#'  \item{\code{raw }}{(sparse) matrix with genes in rows, cells in columns}
#'  \item{\code{totalUMI }}{the colSums of obj$raw}
#'  \item{\code{embed }}{two-dimensional embedding of the cells, provided as data.frame
#'  or tibble with two columns and one row per cell.}
#'  \item{\code{neighbors }}{index matrix with one row per cell and k columns, where
#'  k is the number of nearest neighbors (we recommend 15<k<100, e.g. k=50). 
#'  Here are two ways to get the neighbors index matrix:
#'  \itemize{
#'    \item Use \code{find_knn(featureMatrix)$idx}, where featureMatrix could be
#'    principal components, latent variables or normalized genes (features in
#'     rows, cells in columns).
#'    \item use \code{as(seurat@graphs[["RNA_nn"]], "dgCMatrix")> .1} to extract
#'    the kNN
#'    graph computed on RNA. The \code{> .1} ensures this also works with RNA_snn,
#'     wknn/wsnn or any other
#'    available graph – check with \code{names(seurat@graphs)}. 
#'    }
#'  }
#'  }
