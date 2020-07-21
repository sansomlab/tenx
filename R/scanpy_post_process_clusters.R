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
    make_option(c("--clusters"), default="scanpy.clusters.tsv.gz",
                help="the scanpy cluster assignments"),
    make_option(c("--predefined"), default="none",
                help="a file containing a set of predefined clusters"),
    make_option(c("--algorithm"), default="leiden",
                help="name of the clustering alogorithm"),
    make_option(c("--mincells"), type="integer", default=10,
                help="clusters with fewer cells are set to NA"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

s <- loadSeurat(path=opt$seuratobject)

message("seurat_cluster.R running with default assay: ", DefaultAssay(s))

if(opt$predefined=="none")
    {
        clusters <- read.table(opt$clusters, sep="\t", header=T, as.is=T)
        cluster_ids <- clusters[[opt$algorithm]]

        if(!all(names(cluster_ids)==Cells(s)))
        {    stop("barcodes from scanpy do not match cells from seurat") }

        print(table(cluster_ids))
        x <- table(cluster_ids)

        rejected_clusters <- names(x[x<opt$mincells])
        cluster_ids[cluster_ids %in% rejected_clusters] <- "911"


        print(head(cluster_ids))
        print(table(cluster_ids))
    } else {

        clusters <- read.table(opt$predefined, sep="\t", header=T, as.is=T)

        cluster_ids <- clusters$cluster_id

        }

cluster_ids <- factor(as.numeric(cluster_ids))
names(cluster_ids) <- clusters$barcode


unique_cluster_ids <- unique(cluster_ids)
print(unique_cluster_ids)

print(sprintf("saving a unique list of cluster ids"))
write.table(unique_cluster_ids, file=file.path(opt$outdir,"cluster_ids.tsv"),
            quote=FALSE, col.names = FALSE, row.names = FALSE)

cluster_colors <- gg_color_hue(length(unique_cluster_ids))
print(sprintf("saving the cluster colors (ggplot)"))
write.table(cluster_colors, file=file.path(opt$outdir,"cluster_colors.tsv"),
            quote=FALSE, col.names = FALSE, row.names = FALSE)

print(sprintf("saveRDS of cluster ids"))
saveRDS(cluster_ids, file=file.path(opt$outdir,"cluster_ids.rds"))

cluster_assignments <- data.frame(barcode=as.character(names(cluster_ids)),
                                  cluster_id=as.character(cluster_ids))

write.table(cluster_assignments,
            gzfile(file.path(opt$outdir, "cluster_assignments.tsv.gz")),
            sep="\t", col.names=T, row.names=F, quote=F)


print(sprintf("seurat_cluster.R final default assay: %s", DefaultAssay(s)))

print("Completed")
