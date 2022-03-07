## code to prepare `simulated_umis` dataset goes here
library(Matrix)


# Simulate 3 cell types with different numbers of cells:
n1 = 600
n2 = 200
n3 = 1100
# features:
set.seed(42)
raw <- data.frame(
  CD3E = rpois(n1+n2+n3, lambda=c(rep(1,n1), rep(1,n2), rep(.1,n3))),
  FOXP3= rpois(n1+n2+n3, lambda=c(rep(.001,n1), rep(.5,n2), rep(.001,n3))),
  ICOS = rpois(n1+n2+n3, lambda=c(rep(3,n1), rep(c(3, .1),n2/2), rep(.1,n3))),
  MS4A1= rpois(n1+n2+n3, lambda=c(rep(.1,n1), rep(.1,n2), rep(1,n3))),
  # gene = rpois(n1+n2+n3, lambda=c(rep(1,n1), rep(1,n2), rep(1,n3)))
  totalUMI = round(rlnorm(n1+n2+n3, log(1800), .25))
)
raw <- as(as.matrix(raw), "dgCMatrix")
# with few genes, uwot returns strange embeddings, so I simulate by hand:
set.seed(42)
ump <- data.frame(
  u1 = rnorm(n1+n2+n3, mean = c(rep(-1,n1), rep(c(.5,1.5),n2/2), rep(8,n3)), sd=.5),
  u2 = rnorm(n1+n2+n3, mean = c(rep(5,n1), rep(3.5,n2), rep(1,n3)), sd=.25),
  celltype= c(rep("Th",n1), rep("Treg",n2), rep("B",n3))
)

ump %>% ggplot(aes(u1,u2, col=celltype))+geom_point()+ coord_fixed()
scUtils::feat(ump[,1:2], raw[, "FOXP3"]/raw[,"totalUMI"])
scUtils::feat(ump[,1:2], raw[, "ICOS"]/raw[,"totalUMI"])

#nn <- scpr::get.knn.annoy(as.matrix(sqrt(raw[,-ncol(raw)] / raw[, "totalUMI"])))
tmp <- uwot::umap(as.matrix(sqrt(raw[,-ncol(raw)] / raw[, "totalUMI"])),
                  n_neighbors = 50,
                  ret_nn=TRUE)

simulated_umis <- list(
  # March 2022: I start following genes-in-row conventions that's prevalent in Bioinformatics.
  raw=t(raw),
  neighbors=tmp$nn$euclidean$idx,
  embed    = ump[,1:2],
  celltype = ump$celltype
)

usethis::use_data(simulated_umis, overwrite = TRUE)

# gives weird embedding, I'm too lazy to tweak UMAP so I don't pursue this further:
# ump <- uwot::umap(sqrt(raw[,-ncol(raw)] / raw[, "totalUMI"]))