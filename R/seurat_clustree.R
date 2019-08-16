## Build a clustree showing the relationship between clusterings.

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(tenxutils),
  require(clustree),
  require(xtable)
)

# Options ----

option_list <- list(
    make_option(c("--resolutions"), default=NULL,
                help="a comma separated list of the cluster resolutions"),
    make_option(c("--clusteridfiles"), default=NULL,
                help="a comma separated list of the rds files with the cluster assignments"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")

    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

resolutions = unlist(strsplit(opt$resolutions,","))
cluster_id_files = unlist(strsplit(opt$clusteridfiles,","))
names(cluster_id_files) <- resolutions

# construct the index of cluster membership
begin = T
for(resolution in resolutions)
{
    tmp <- readRDS(cluster_id_files[resolution])

    if(begin==T) {
        clust_index = data.frame(row.names=names(tmp))
        begin = F
    }
    clust_index[[resolution]] = tmp[rownames(clust_index)]
}

colnames(clust_index) <- paste0("R",colnames(clust_index))

# s@meta.data <- clust_index

gp <- clustree(clust_index,"R")

save_ggplots(
    file.path(opt$outdir, "clustree"),
    gp=gp,
    width=6, height=8
    )

message("Completed")
