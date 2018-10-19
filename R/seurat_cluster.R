## Cluster the single cells in a given Seurat object

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
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--components"), type="integer", default=10,
                help="number of principle components to use"),
    make_option(c("--resolution"), type="double", default=1,
                help="cluster resolution"),
    make_option(c("--algorithm"), type="integer", default=3,
                help="1=original Louvain, 2=Louvain multilevel, 3=SLM"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use in construction of SNN graph. (e.g. 'pca', 'ica')")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# set the run specs
run_specs <- paste(opt$components,opt$resolution,opt$algorithm,opt$testuse,sep="_")


message(sprintf("readRDS: %s", opt$seuratobject))
s <- readRDS(opt$seuratobject)

## The FindClusters() function implements the procedure, and contains a
## resolution parameter that sets the 'granularity' of the downstream clustering,
## with increased values leading to a greater number of clusters.
## We find that setting this parameter between 0.6-1.2 typically returns
## good results for single cell datasets of around 3K cells. Optimal resolution
## often increases for larger datasets. The clusters are saved in the object\@ident slot.

message(sprintf("FindClusters"))
s <- FindClusters(s,
                  reduction.type=opt$reductiontype,
                  dims.use = 1:opt$components,
                  resolution = opt$resolution,
                  algorithm = opt$algorithm,
                  print.output = 0,
                  save.SNN = F)

cluster_ids <- s@ident

nclusters <- length(unique(cluster_ids))

message(sprintf("write"))
write(nclusters, file=file.path(opt$outdir,"nclusters.txt"))
message(sprintf("saveRDS"))
saveRDS(cluster_ids, file=file.path(opt$outdir,"cluster_ids.rds"))

## Visualise the relationship between the clusters.

draw_tree <- function() { BuildClusterTree(s) }

save_plots(
    file.path(opt$outdir, "cluster.dendrogram"),
    plot_fn=draw_tree,
    width=8, height=5
    )




message("Completed")
