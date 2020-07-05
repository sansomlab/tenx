## Cluster the single cells in a given Seurat object

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(tenxutils),
  require(xtable)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--clusterids"), default="cluster_ids.rds",
                help="A seurat object after PCA"),
    make_option(c("--usesigcomponents"), default=FALSE,
                help="use significant principle component"),
    make_option(c("--components"), type="integer", default=10,
                help="if usesigcomponents is FALSE, the number of principle components to use"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use in construction of SNN graph. (e.g. 'pca', 'ica')")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

message(sprintf("readRDS: %s", opt$seuratobject))
s <- readRDS(opt$seuratobject)

message("seurat_cluster.R running with default assay: ", DefaultAssay(s))

Idents(s) <- readRDS(opt$clusterids)

print(Idents(s))

## Get the components
if(opt$usesigcomponents)
{
    comps <- getSigPC(s)
    message("using the following pcas:")
    print(comps)
} else {
    comps <- 1:as.numeric(opt$components)
}


## Visualise the relationship between the clusters.

## 1. Using the Seurat default function
## which computes distance between clusters averages
## using the expression levels of the variable genes

if(DefaultAssay(s) != "SCT" && length(VariableFeatures(s)) > 0)
{
    ## see: https://github.com/satijalab/seurat/issues/1677

    draw_tree <- function() { plot(Tool(BuildClusterTree(s),slot='BuildClusterTree')) }

} else {
    ## draw an empty plot with an error message

    draw_tree <- function()
    {
    plot.new()
    text(0.5,0.5,"BuildClusterTree is not compatible with sctransform")
    }
}

save_plots(
    file.path(opt$outdir, "cluster.dendrogram"),
    plot_fn=draw_tree,
    width=8, height=5)


## 2. By the median pair-wise pearson correlation
## of cells in the clusters.
cluster_cor <- clusterCor(s,
                          dr_type=opt$reductiontype,
                          comps,
                          cor_method="pearson",
                          cluster_average=FALSE)

den <- as.dendrogram(hclust(as.dist(1 - cluster_cor)))

draw_cor_tree <- function() { plot(den) }

save_plots(
    file.path(opt$outdir, "cluster.correlation.dendrogram"),
    plot_fn=draw_cor_tree,
    width=8, height=5
    )

## The diagonal of the correlation matrix is interesting
## because it is a measure of with-cluster heterogenity
## so we write the table out
print(
    xtable(sprintfResults(as.data.frame(cluster_cor)), caption="Pairwise correlations between clusters"),
    file=file.path(opt$outdir, "cluster.pairwise.correlations.tex")
    )

message("seurat_cluster.R final default assay: ", DefaultAssay(s))

message("Completed")
