## Compute reduced tSNE dimensions and save to a table
## along with the cell metadata

# Libraries ----

stopifnot(
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(optparse),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object containing reduced dimensions"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--annotation"), default="none",
                help="A file containing the mapping of gene_id, gene_name and seurat_id"),
    make_option(c("--usesigcomponents"), default=FALSE,
                help="use significant principle component"),
    make_option(c("--components"), type="integer", default=10,
                help="number of reduced dimension components to use"),
    make_option(c("--perplexity"), type="integer", default=30,
                help="the value of the perplexity hyper-parameter"),
    make_option(c("--maxiter"), type="integer", default=5000,
                help="maximum number of iterations"),
    make_option(c("--fast"), default=TRUE,
                help="if true the faster but less flexible Barnes-hut implementation is used"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use in construction of SNN graph. (e.g. 'pca', 'ica', 'zinbwave', etc)"),
    make_option(c("--outfile"), default="tsne.txt",
                help="the file to which the tSNE coordinates will be written")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


message("readRDS")
s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
Idents(s) <- cluster_ids
print(head(Idents(s)))

message("seurat_tsne.R running with default assay: ", DefaultAssay(s))

## check that the perplexity value specified is sensible
ncells <- ncol(GetAssayData(object = s))
message("no. cells:", ncells)
message("perplexity:", opt$perplexity)

## get the principle components to use
if(opt$usesigcomponents)
{
    comps <- getSigPC(s)
    message("using the following pcas:")
    print(comps)
} else {
    comps <- 1:as.numeric(opt$components)
}

if(opt$perplexity > floor(ncells/5))
{
    message("Perplexity > floor(ncells/5), skipping")
} else {

    ## run the tSNE analysis
    message("RunTSNE")
    s <- RunTSNE(s,
                 reduction = opt$reductiontype,
                 dims = comps,
                 perplexity = opt$perplexity,
                 max_iter = opt$maxiter,
                 do.fast = as.logical(opt$fast))

    ## extract the tSNE coordinates from the seurat object
    tsne <- as.data.frame(s@reductions$tsne@cell.embeddings)
    tsne$cluster <- Idents(s)[rownames(tsne)]

    plot_data <- merge(tsne, s[[]], by=0)

    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL

    plot_data$barcode <- row.names(plot_data)

    ## save the annotated tSNE data frame
    write.table(plot_data, opt$outfile,
                sep="\t", quote=FALSE, row.names=FALSE)

message("Completed")

    }

message("seurat_tsne.R final default assay: ", DefaultAssay(s))
