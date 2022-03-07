# Documents the data set generated in data-raw/simulated_umis.R
# (see source code of cellpypes package).



#' Simulated scRNAseq data
#'
#' This data serves to develop cellpypes and to illustrate its functionality.
#' I made it up entirely.
#'
#' @format A list with 4 entries:
#' \describe{
#'   \item{raw}{Raw (unnormalized) UMI counts for a handful of genes, last row are totalUMI. }
#'   \item{neighbors}{Indices of each cell's 50 nearest neighbors.}
#'   \item{embed}{Simulated UMAP embedding.}
#'   \item{celltype}{Cell type label that I used to simulate the data.}
#' }
#' @source Very simple simulation (c.f. data-raw/simulated_umis.R in source code).
"simulated_umis"
