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
    make_option(c("--seuratgraphs"), default="begin.Robj",
                help="An rds object containing the seurat nn and snn graphs"),

    ## make_option(c("--usesigcomponents"), default=FALSE,
    ##             help="use significant principle component"),
    ## make_option(c("--components"), type="integer", default=10,
    ##             help="if usesigcomponents is FALSE, the number of principle components to use"),
    make_option(c("--predefined"), default=NULL,
                help="An rds file containing a named vector of predefined cell cluster assignments"),
    make_option(c("--resolution"), type="double", default=1,
                help="cluster resolution"),
    make_option(c("--algorithm"), type="integer", default=3,
                help="1=original Louvain, 2=Louvain multilevel, 3=SLM"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

message(sprintf("readRDS: %s", opt$seuratobject))
s <- readRDS(opt$seuratobject)
s@graphs <- readRDS(opt$seuratgraphs)

message("seurat_cluster.R running with default assay: ", DefaultAssay(s))

## The FindClusters() function implements the procedure, and contains a
## resolution parameter that sets the 'granularity' of the downstream clustering,
## with increased values leading to a greater number of clusters.
## We find that setting this parameter between 0.6-1.2 typically returns
## good results for single cell datasets of around 3K cells. Optimal resolution
## often increases for larger datasets. The clusters are saved in the object\@ident slot.

## if(opt$usesigcomponents)
## {
##     comps <- getSigPC(s)
##     message("using the following pcas:")
##     print(comps)
## } else {
##     comps <- 1:as.numeric(opt$components)
## }

## if(toupper(opt$reductiontype)=="PCA")
## { if(nrow(s@reductions$pca@jackstraw@overall.p.values) > 0)
##   {
##     ## Make a table of the retained principle components
##     x <- as.data.frame(s@reductions$pca@jackstraw@overall.p.values)
##     x$p.adj <- p.adjust(x$Score, method="BH")
##     x$significant <- "no"
##     x$significant[x$p.adj < 0.05] <- "yes"
##     x <- x[x$PC %in% comps,]
##     x$sdev <- s@reductions$pca@stdev[x$PC]

##     print(
##         xtable(sprintfResults(x), caption=paste("Table of the selected (n=",
##                                                 nrow(x),
##                                                 ") principle components",
##                                                 sep="")),
##         file=file.path(opt$outdir, "selected.principal.components.tex")
##     )
##   }
## }


if(is.null(opt$predefined))
{


    message(sprintf("FindClusters"))


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



unique_cluster_ids <- sort(unique(cluster_ids))

message(sprintf("saving a unique list ofe cluster ids"))
write.table(unique_cluster_ids, file=file.path(opt$outdir,"cluster_ids.txt"),
            quote=FALSE, col.names = FALSE, row.names = FALSE)

cluster_colors <- gg_color_hue(length(unique_cluster_ids))
message(sprintf("saving the cluster colors (ggplot)"))
write.table(cluster_colors, file=file.path(opt$outdir,"cluster_colors.txt"),
            quote=FALSE, col.names = FALSE, row.names = FALSE)

message(sprintf("saveRDS"))
saveRDS(cluster_ids, file=file.path(opt$outdir,"cluster_ids.rds"))

cluster_assignments <- data.frame(barcode=as.character(names(cluster_ids)),
                                  cluster_id=as.character(cluster_ids))

write.table(cluster_assignments,
            gzfile(file.path(opt$outdir, "cluster_assignments.txt.gz")),
            sep="\t", col.names=T, row.names=F, quote=F)


message("seurat_cluster.R final default assay: ", DefaultAssay(s))

message("Completed")
