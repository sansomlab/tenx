## Run SingleR in a given Seurat object

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(dplyr),
  require(tenxutils),
  require(SingleR),
  require(scater),
  require(tidyr),
  require(ggforce),
  require(MASS),
  require(viridis)
)

# Options ----

option_list <- list(
  make_option(c("--seuratobject"), default="begin.Robj",
              help="A seurat object after PCA"),
  make_option(c("--reference"),
              help="reference dataset for celltype assignment; see https://bioconductor.org/packages/3.11/bioc/vignettes/SingleR/inst/doc/SingleR.html#5_available_references" ),
  make_option(c("--workers"), default=10,
              help="number of parallel processes to use"),
  make_option(c("--show_annotation_in_plots"), default=NULL,
              help="Column names from the metadata slot to show in the annotation of the output plots"),
  make_option(c("--outdir"), default="seurat.out.dir",
              help="outdir")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

message(sprintf("readRDS: %s", opt$seuratobject))
s <- readRDS(opt$seuratobject)
a <- ifelse("SCT" %in% names(s), yes = "SCT", no = "RNA")

sce <- as.SingleCellExperiment(s, assay = a) # Seems to be using SCT assay (so does it take the defualt?)

cat("Retrieving SummarizedExperiment object for reference datset:", opt$reference, "\n\n")
ref.se <- match.fun(opt$reference)
ref.se <- ref.se()

# Run on common genes
common <- intersect(rownames(sce), rownames(ref.se))
ref.se <- ref.se[common,]
sce <- sce[common,]
cat("\n\nUsing", format(length(common), big.mark = ","), "genes for prediction \n\n")

# Predict cell identity
cat("Predicting cell identity \n")

multicoreParam <- MulticoreParam(workers = opt$workers)

pred <- SingleR(test = sce,
                ref = ref.se,
                labels = ref.se$label.main,
                BPPARAM = multicoreParam)

saveRDS(pred,
        file.path(opt$outdir,"predictions.rds"))

labels <- data.frame(barcode=rownames(pred),
                     first.labels=pred$first.labels,
                     labels=pred$labels,
                     pruned.labels=pred$pruned.labels)

write.table(labels,
            gzfile(file.path(opt$outdir,"labels.tsv.gz")),
            col.names=TRUE, row.names=FALSE,
            sep="\t", quote=FALSE)


message("Completed\n\n")
sessionInfo()
