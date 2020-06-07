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
#  make_option(c("--reference"),
#              help="reference dataset for celltype assignment; see https://bioconductor.org/packages/3.11/bioc/vignettes/SingleR/inst/doc/SingleR.html#5_available_references" ),
  make_option(c("--predictions"), default = NULL,
              help="rds object containig the result of the call to singleR"),
  make_option(c("--show_annotation_in_plots"), default=NULL,
              help="Column names from the metadata slot to show in the annotation of the output plots"),
  make_option(c("--outdir"), default="seurat.out.dir",
              help="outdir"),
  make_option(c("--plotformat"), default="png",
              help="pdf, png, or both")

)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

message(sprintf("readRDS: %s", opt$seuratobject))
s <- readRDS(opt$seuratobject)
a <- ifelse("SCT" %in% names(s), yes = "SCT", no = "RNA")

sce <- as.SingleCellExperiment(s, assay = a) # Seems to be using SCT assay (so does it take the defualt?)

# read in the predictions
pred <- readRDS(opt$predictions)

common <- intersect(rownames(sce), rownames(pred))
sce <- sce[common,]

# Plot label assignment scores
cat("Plotting prediction results \n")
annotation_frame <- as.data.frame(colData(sce))
qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mito")

if (!is.null(opt$show_annotation_in_plots)) {
    meta_cols <- unlist(strsplit(opt$show_annotation_in_plots, split = ","))

  keep_annotation <- c(qc_cols, meta_cols)
} else {
  keep_annotation <- qc_cols
}

keep_annotation <- keep_annotation[keep_annotation %in% colnames(annotation_frame)]

annotation_frame <- annotation_frame[,keep_annotation]

do_plot <- function() {
plotScoreHeatmap(pred, show.labels = TRUE,
                 annotation_col=annotation_frame)
}

save_plots(paste0(opt$outdir, "/singleR_score_heatmap"),
           plot_fn=do_plot,
           width = 12,
           height = 8)


message("Completed\n\n")
