## Title ----
##
## Aggregate UMI counts across cells within cluster to form pseudobulks.
##
## Description ----
##
## A script to that use cluster information to aggregate UMI counts
## of cells within each cluster, and stores the resulting matrix of
## pseudo-bulks UMI counts to an RDS file.
##
## Usage ----
##
## statement = '''
##                &> %(outfile)s
##             '''

message("aggregate_cluster.R")
timestamp()

# Libraries ----

stopifnot(
    require(optparse),
    require(Seurat),
    require(dplyr),
    require(Matrix),
    require(reshape2),
    require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
        c("--tenxdir"),
        help="Path to the input 10x matrix directory."),
    make_option(
        c("--clusterids"),
        help="Path to an RDS file containing a named vector of cluster assignments."),
    make_option(
        c("--outfile"),
        default="cluster_counts.rds",
        help="RDS file to write the aggregated count matrix.")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Input data ----

## Matrix
# TODO: used DropletUtils package instead
cat("Importing matrix from: ", opt$tenxdir, " ... ")
matrixUMI <- Read10X(opt$tenxdir)
cat("Done.\n")
cat(
    "Input matrix size:",
    sprintf("%i rows/genes, %i columns/cells\n", nrow(matrixUMI), ncol(matrixUMI))
)

## Clusters
cat("Importing cluster information from:", opt$clusterids, "...")
clusterIds <- readRDS(opt$clusterids)
cat("Done.\n")
cat("Cluster assignments:", length(clusterIds), "\n")

stopifnot(is.factor(clusterIds))
stopifnot(!is.null(names(clusterIds)))
stopifnot(
    all(names(clusterIds) %in% colnames(matrixUMI))
)

# Compute pseudo-bulks ----

clusterCountMatrix <- matrix(
    nrow = nrow(matrixUMI),
    ncol = nlevels(clusterIds),
    dimnames = list(
        gene_ids = rownames(matrixUMI),
        cluster = paste("cluster", levels(clusterIds), sep="_")
    ))
for (clusterId in levels(clusterIds)){
    cat("Processing cluster:", clusterId, "...")
    cellsInCluster <- names(clusterIds)[which(clusterIds == clusterId)]
    pseudobulkCount <- rowSums(matrixUMI[, cellsInCluster])
    clusterCountMatrix[, paste("cluster", clusterId, sep = "_")] <-
        pseudobulkCount[rownames(clusterCountMatrix)]
    cat("Done.\n")
}

# Write out matrix ----
cat("Writing pseudobulk matrix to RDS file ...")
saveRDS(clusterCountMatrix, opt$outfile)
cat("Done.\n")

timestamp()
message("Completed")
