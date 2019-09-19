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
    make_option(c("--usesigcomponents"), default=FALSE,
                help="use significant principle component"),
    make_option(c("--components"), type="integer", default=10,
                help="if usesigcomponents is FALSE, the number of principle components to use"),
    make_option(c("--predefined"), default=NULL,
                help="An rds file containing a named vector of predefined cell cluster assignments"),
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


message("seurat_cluster.R running with default assay: ", DefaultAssay(s))

## The FindClusters() function implements the procedure, and contains a
## resolution parameter that sets the 'granularity' of the downstream clustering,
## with increased values leading to a greater number of clusters.
## We find that setting this parameter between 0.6-1.2 typically returns
## good results for single cell datasets of around 3K cells. Optimal resolution
## often increases for larger datasets. The clusters are saved in the object\@ident slot.

if(opt$usesigcomponents)
{
    comps <- getSigPC(s)
    message("using the following pcas:")
    print(comps)
} else {
    comps <- 1:as.numeric(opt$components)
}

if(toupper(opt$reductiontype)=="PCA")
{
    ## Make a table of the retained principle components
    x <- as.data.frame(s@reductions$pca@jackstraw@overall.p.values)
    x$p.adj <- p.adjust(x$Score, method="BH")
    x$significant <- "no"
    x$significant[x$p.adj < 0.05] <- "yes"
    x <- x[x$PC %in% comps,]
    x$sdev <- s@reductions$pca@stdev[x$PC]

    print(
        xtable(sprintfResults(x), caption=paste("Table of the selected (n=",
                                                nrow(x),
                                                ") principle components",
                                                sep="")),
        file=file.path(opt$outdir, "selected.principal.components.tex")
    )
}


if(is.null(opt$predefined))
{


    message(sprintf("FindClusters"))
    s <- FindNeighbors(s,
                       reduction.type = "pca",
                       dims = comps)


    s <- FindClusters(s,
                      resolution = opt$resolution,
                      algorithm = opt$algorithm)


    cluster_ids <- Idents(s) #

} else {

    user_ids <- readRDS(opt$predefined)

    ## check that the barcodes match.
    if(!all(names(user_ids) %in% Cells(s)) | length(user_ids)!= length(Cells(s)))
    {
        stop(paste("The cluster assignments supplied by the user do not contain ",
                   "the same cells as the seurat object", sep=""))
    }

    cluster_ids <- user_ids[Cells(s)]

    Idents(s) <- cluster_ids

    ## reorder the cells to match the seurat object

}



nclusters <- length(unique(cluster_ids))

message(sprintf("write"))
write(nclusters, file=file.path(opt$outdir,"nclusters.txt"))
message(sprintf("saveRDS"))
saveRDS(cluster_ids, file=file.path(opt$outdir,"cluster_ids.rds"))

## Visualise the relationship between the clusters.

## 1. Using the Seurat default function
## which computes distance between clusters averages
## using the expression levels of the variable genes

if(DefaultAssay(s) != "SCT")
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
