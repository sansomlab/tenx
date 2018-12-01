## Compute reduced UMAP dimensions and save to a table
## along with the cell metadata

# Libraries ----

stopifnot(
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--annotation"), default="none",
                help="A file containing the mapping of gene_id, gene_name and seurat_id"),
    make_option(c("--components"), type="integer", default=10,
                help="number of principle components to use"),
    make_option(c("--nneighbors"), type="integer", default=30L,
                help="number of neighboring points used for local approximation of the manifold structure"),
    make_option(c("--mindist"), type="double", default=0.3,
                help="controls how tightly the embedding is allowed to compress points together"),
    make_option(c("--metric"), type="character", default="correlation",
                help="the choice of metric used to measure distance in the input space"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction slot to use as the input to UMAP (e.g. 'pca', 'ica')"),
    make_option(c("--outfile"), default="umap.txt",
                help="the file to which the UMAP coordinates will be written")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


message("readRDS")
s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
s@ident <- cluster_ids

## run UMAP
message("RunUMAP")
s <- RunUMAP(s,
             cells.use = NULL,
             reduction.use = opt$reductiontype,
             dims.use = 1:opt$components,
             genes.use = NULL,
             assay.use = "RNA",
             max.dim = 2L,
             reduction.name = "umap",
             reduction.key = "UMAP",
             n_neighbors = opt$nneighbors,
             min_dist = opt$mindist,
             metric = opt$metric,
             seed.use = 42)

## extract the UMAP coordinates from the seurat object
umap <- as.data.frame(s@dr$umap@cell.embeddings)
umap$cluster <- s@ident

plot_data <- merge(umap, s@meta.data, by=0)

rownames(plot_data) <- plot_data$Row.names
plot_data$Row.names <- NULL

plot_data$barcode <- row.names(plot_data)

## save the annotated tSNE data frame
write.table(plot_data, opt$outfile,
            sep="\t", quote=FALSE, row.names=FALSE)

message("Completed")
