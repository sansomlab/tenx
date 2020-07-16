## Title ----
##
## Visualise gene expression levels on reduced dimensions plots
##
## Description ----
##
## This script visualises per-cell gene expression levels on a 2D
## reduced dimensions plot.
##
## Details ----
##
## Expression levels are taken from the "data" slot,
## and truncated at the 95th expression percentile
##
## Usage ----
##
## See options.

# Libraries ----

stopifnot(
  require(optparse),
  require(ggplot2),
  require(reshape2),
  require(Seurat),
  require(tenxutils),
  require(BiocParallel),
  require(Matrix)
)

# Options ----

option_list <- list(
    make_option(
      c("--features"),
      default="none",
      help="A tab-delimited text file containing gene_id (or gene if s@misc$gene_id is not set) and gene_name"
      ),
    make_option(
      c("--seuratobject"),
      default="none",
      help="The seurat object (e.g. begin.rds)"
    ),
    make_option(
      c("--seuratassay"),
      default="RNA",
      help="The seurat assay to pull the expression data from"
    ),
    make_option(
      c("--slot"),
      default="None",
      help="The name of the slot to pull the data from"
      ),
    make_option(
      c("--outfile"),
      default="assay.data.tsv.gz",
      help="the file to write the data to")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


if (endsWith(opt$seuratobject, ".rds")) {
  message(sprintf("readRDS: %s", opt$seuratobject))
  s <- readRDS(opt$seuratobject)
} else {
  message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
  stopifnot(require(SeuratDisk))
  s <- LoadH5Seurat(opt$seuratobject)
}

## set the default assay
message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("plot_rdims_gene.R running with default assay: ", DefaultAssay(s))
data <- GetAssayData(object = s, slot = opt$slot)  #data!

# remove the Seurat object
rm(s)

## read in the table containing the genes to visualise
features <- read.table(opt$features,
                       header=TRUE,
                       sep="\t",
                       as.is=TRUE)

## subset the gene expression matrix
exprs <- as.matrix(data[features$gene,, drop=F])

saveRDS(exprs, opt$outfile)
#write.table(exprs, gzfile(opt$outfile), sep="\t", col.names=T, quote=F)

message("completed")
