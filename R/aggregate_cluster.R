## Cluster the single cells in a given Seurat object

library(optparse)

## deal with the options
option_list <- list(
    make_option(c("--tenxdir"), default="must_specify",
                help="location of 10x matrix"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--outfile"), default="cluster_counts.rds",
                help="RDS file that stores the aggregated count matrix")
    )

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with options:")
print(opt)

# load required libraries
stopifnot(
    require(Seurat),
    require(dplyr),
    require(Matrix),
    require(reshape2)
    require(tenxutils)
)

# import 10x raw counts
message("Reading 10x dir")
data10x <- Read10X(opt$tenxdir)
message("10x matrix size:")
print(dim(data10x))

# import cluster assignments
clusterIds <- readRDS(opt$clusterids)

stopifnot(
    all(names(clusterIds) %in% colnames(data10x))
)

clusterCountMatrix <- matrix(
  nrow = nrow(data10x),
  ncol = nlevels(clusterIds),
  dimnames = list(
    gene_ids = rownames(data10x),
    cluster = paste("cluster", levels(clusterIds), sep="_")
))
for (clusterId in levels(clusterIds)){
  message("cluster: ", clusterId)
  cellsInCluster <- names(clusterIds)[which(clusterIds == clusterId)]
  pseudobulkCount <- rowSums(data10x[,cellsInCluster])
  clusterCountMatrix[,paste("cluster", clusterId, sep = "_")] <- pseudobulkCount[rownames(clusterCountMatrix)]
}

saveRDS(clusterCountMatrix, opt$outfile)
message("Completed")
