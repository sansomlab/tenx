#' ---
#' title: "Distribution of UMIs, number genes and fraction of mitochondrial UMIs from all-datasets.dir/all/"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/lkoneva/tenx_qctest/Rmd/qc_scatterplots.yml"
#'  fig_path: "fig.dir/"
#'  log_filename: "qc_scatterplots.log"
#'  sample: NA
#' ---
#' ---
#' #'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------

stopifnot(require(optparse))
stopifnot(require(yaml))
stopifnot(require(ggplot2))
stopifnot(require(dplyr))
stopifnot(require(RSQLite))
stopifnot(require(Matrix))
stopifnot(require(reshape2))
stopifnot(require(future))
stopifnot(require(grid))
stopifnot(require(tenxutils))
stopifnot(require(viridis))
stopifnot(require(gridExtra))
stopifnot(require(ggExtra))
stopifnot(require(knitr))
stopifnot(require(rtracklayer))
stopifnot(require(DropletUtils))
stopifnot(require(futile.logger))

# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               include = FALSE,
               fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:
default_options <- list(
  # Location of the code 
  "tenx_dir" = "",
  
  # Name of folders to output rds, pdf and png files
  "outdir" = "",
  
  # Path to cellranger all-datasets.dir/all/ 
  "matrixpath" = "",

  # path to gtf file
  "transcriptome" = "cellranger_transcriptome",

  # current sample
  "sample" = "", 
  
  # parameters for filtering in pipeline_seurat
  # shown as blue dot lines on scatterplots
  "qc_initial_mingenes" = 200,
  "qc_mingenes" = 500,
  "qc_minpercentmito" = 0,
  "qc_maxpercentmito" = 0.05
)

# here the yaml can also be read in from the default location in code 
# directory
options <- read_yaml(params$task_yml)
#options <- read_yaml("/gfs/work/lkoneva/07_monaco_plaque/analysis/aneurysm/cellranger3_qctest/qc.dir/qc_scatterplots.yml")

# Update the default options
if(!is.null(options)) {
 opt <- utils::modifyList(default_options, options)
} else{
 opt <- default_options
}

# make output directories per sample as in pipeline_cellranger.py
opt[["sample"]] = params$sample
out = file.path(opt$outdir, paste(opt$sample, "_plot.dir", sep = ""))
opt[["outdir"]] = out

# Prepare the logger and the log file location
# the logger is by default writing to ROOT logger and the INFO threshold level
flog.threshold(INFO)
# now set to append mode
flog.appender(appender.file(params$log_filename))

flog.info("Running with options:", opt, capture = TRUE)
flog.info("\n")

## ######################################################################### ##
## ###################### (0) Read in data ################################# ##
## ######################################################################### ##

flog.info("Reading gtf file ...")
gtfpath <- file.path(opt$transcriptome, "genes/genes.gtf")
gtfGeneInfo <- import.gff(gtfpath, feature.type="gene")
names(gtfGeneInfo) <- gtfGeneInfo$gene_id
flog.info("Done.\n")

flog.info("Reading metadata ...")
metadatapath <- file.path(opt$matrixpath, "metadata.tsv.gz")
metadata <- read.table(gzfile(metadatapath), sep="\t", header=TRUE, as.is=TRUE)
flog.info("Done.\n")

# function to read matrices 
qc_mat <- function(matrixpath){
  sce10x <- read10xCounts(matrixpath, type = "auto", col.names=TRUE)
  rowRanges(sce10x) <- gtfGeneInfo[rownames(sce10x)]
  sce10x$sample <- sce10x$Sample
  sce10x$umis <- colSums(assay(sce10x, "counts"))
  sce10x$umis_mt <- colSums(assay(sce10x[seqnames(sce10x) == "MT", ], "counts"))
  sce10x$genes <- colSums(assay(sce10x, "counts") > 0)
  tableExport <- colData(sce10x)[, c("sample", "umis", "umis_mt", "genes")]
}

flog.info("Creating QC table ...")
tableData <- as.data.frame(qc_mat(opt$matrixpath))
tableData$sample <- metadata$sample_id

## fraction of mito umis
flog.info("Preprocessing data for ggplot ...")
tableData$umis_mt_percent <- tableData$umis_mt / tableData$umis
flog.info("Done.\n")

# select data for sample N
tmp <- tableData %>% filter(tableData$sample == opt$sample)

## ######################################################################### ##
## ###################### (i) Genes vs UMIs ################################ ##
## ######################################################################### ##

flog.info("Creating scatter plot...")

#' QC plots for sample: `r opt$sample`
#' 

p1 <- ggplot(tmp, aes(x=genes, y=umis, color=umis_mt_percent)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_bw() +
  xlab("Number of genes") +
  ylab("Total UMIs") +
  # add thresholds
  geom_vline(xintercept=c(100,200,500,1000), linetype="dashed", color = "red", size=0.1) +
  geom_hline(yintercept=c(200,500,1000,5000), linetype="dashed", color = "red", size=0.1) +  
  geom_vline(xintercept=c(opt$qc_initial_mingenes, opt$qc_mingenes), linetype="dashed", color = "blue", size=0.2) +
  # legends
  labs(color = "Fraction mitochondrial counts") +
  theme(legend.position = "bottom",  legend.key.width = unit(1.5,"cm"), legend.key.height = unit(0.2,"cm")) +
  guides(color = guide_colorbar(title.position = "top"))
p1 <- ggMarginal(p1, type="density")

save_ggplots(file.path(opt$outdir, "genes_vs_umis"), p1)

#' Number Genes vs Total UMIs colored by % Mito.
#+ plot_genes_umis, include=TRUE, fig.height=5
print(p1)
#+ include=FALSE
#'

## ######################################################################### ##
## ################### (ii) log10(Genes) vs log10(UMIs) #################### ##
## ######################################################################### ##

flog.info("Creating scatter plot...")

p2 <- ggplot(tmp, aes(x=genes, y=umis, color=umis_mt_percent)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_bw() +
  xlab("log10(Number of genes)") +
  ylab("log10(Total UMIs)") +
  # transform data
  scale_x_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
  scale_y_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
  # add thresholds
  geom_vline(xintercept=c(100,200,500,1000), linetype="dashed", color = "red", size=0.1) +
  geom_hline(yintercept=c(200,500,1000,5000), linetype="dashed", color = "red", size=0.1) +
  geom_vline(xintercept=c(opt$qc_initial_mingenes, opt$qc_mingenes), linetype="dashed", color = "blue", size=0.2) +
  # legends
  labs(color = "Fraction mitochondrial counts") +
  theme(legend.position = "bottom",  legend.key.width = unit(1.5,"cm"), legend.key.height = unit(0.2,"cm")) +
  guides(color = guide_colorbar(title.position = "top"))
p2 <- ggMarginal(p2, type="density")

save_ggplots(file.path(opt$outdir, "log10genes_vs_log10umis"), p2)

#' log10(Number Genes) vs log10(Total UMIs) colored by % Mito.
#+ plot_log10genes_log10umis, include=TRUE, fig.height=5
print(p2)
#+ include=FALSE
#'


## ######################################################################### ##
## ####################### (iii) Number Genes vs % Mito #################### ##
## ######################################################################### ##

flog.info("Creating scatter plot...")

p3 <- ggplot(tmp, aes(x=genes, y=umis_mt_percent, color=umis)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_bw() +
  xlab("Number of genes") +
  ylab("Fraction of mitochondrial UMIs") +
  # add thresholds
  geom_hline(yintercept=c(0.2,0.5), linetype="dashed", color = "red", size=0.2) +
  geom_vline(xintercept=c(100,200,500,1000), linetype="dashed", color = "red", size=0.2) +
  geom_hline(yintercept=c(opt$qc_minpercentmito, opt$qc_maxpercentmito), linetype="dashed", color = "blue", size=0.2) +
  geom_vline(xintercept=c(opt$qc_initial_mingenes, opt$qc_mingenes), linetype="dashed", color = "blue", size=0.2) +
  # legend
  labs(color = "Count depth (number of UMIs)") +
  theme(legend.position = "bottom",  legend.key.width = unit(1.5,"cm"), legend.key.height = unit(0.2,"cm")) +
  guides(color = guide_colorbar(title.position = "top"))
p3 <- ggMarginal(p3, type="density")

save_ggplots(file.path(opt$outdir, "10genes_vs_Mito"), p3)

#' Number Genes vs Fraction of mitochondrial UMIs colored by UMIs
#+ plot_genes_Mito, include=TRUE, fig.height=5
print(p3)
#+ include=FALSE
#'

## ######################################################################### ##
## ################# (iv) log10(Genes) vs sqrt(% Mito) ##################### ##
## ######################################################################### ##

flog.info("Creating scatter plot...")

p4 <- ggplot(tmp, aes(x=genes, y=umis_mt_percent, color=umis)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_bw() +
  xlab("log10(Number of genes)") +
  ylab("sqrt(Fraction of mitochondrial UMIs)") +
  # transform data
  scale_x_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
  scale_y_sqrt(breaks=c(0.00,0.05,0.20,0.25,0.50,0.75,1.00), labels=c(0.00,0.05,0.20,0.25,0.50,0.75,1.00)) +
  # add thresholds
  geom_hline(yintercept=c(0.2,0.5), linetype="dashed", color = "red", size=0.2) +
  geom_vline(xintercept=c(100,200,500,1000), linetype="dashed", color = "red", size=0.2) +
  geom_hline(yintercept=c(opt$qc_minpercentmito, opt$qc_maxpercentmito), linetype="dashed", color = "blue", size=0.2) +
  geom_vline(xintercept=c(opt$qc_initial_mingenes, opt$qc_mingenes), linetype="dashed", color = "blue", size=0.2) +
  # legend
  labs(color = "Count depth (number of UMIs)") +
  theme(legend.position = "bottom",  legend.key.width = unit(1.5,"cm"), legend.key.height = unit(0.2,"cm")) +
  guides(color = guide_colorbar(title.position = "top"))
p4 <- ggMarginal(p4, type="density")

save_ggplots(file.path(opt$outdir, "log10genes_vs_sqrtMito"), p4)

#' log10(Number Genes) vs sqrt(Fraction of mitochondrial UMIs) colored by UMIs
#+ plot_log10genes_sqrtMito, include=TRUE, fig.height=5
print(p4)
#+ include=FALSE
#'

## ######################################################################### ##
## ###################### (v) Total UMIs vs % Mito ######################### ##
## ######################################################################### ##

flog.info("Creating scatter plot...")

p5 <- ggplot(tmp, aes(x=umis, y=umis_mt_percent, color=genes)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_bw() +
  xlab("Total UMIs") +
  ylab("Fraction of mitochondrial UMIs") +
  # add thresholds
  geom_hline(yintercept=c(0.2,0.5), linetype="dashed", color = "red", size=0.2) +
  geom_vline(xintercept=c(200,500,1000,5000), linetype="dashed", color = "red", size=0.2) +
  geom_hline(yintercept=c(opt$qc_minpercentmito, opt$qc_maxpercentmito), linetype="dashed", color = "blue", size=0.2) +
  # legend
  labs(color = "Number of genes") +
  theme(legend.position = "bottom",  legend.key.width = unit(1.5,"cm"), legend.key.height = unit(0.2,"cm")) +
  guides(color = guide_colorbar(title.position = "top"))
p5 <- ggMarginal(p5, type="density")

save_ggplots(file.path(opt$outdir, "log10umis_vs_sqrtMito"), p5)

#' Total UMIs vs Fraction of mitochondrial UMIs colored by Number Genes
#+ plot_umis_Mito, include=TRUE, fig.height=5
print(p5)
#+ include=FALSE
#'

## ######################################################################### ##
## #################### (vi) log10(UMIs) vs % sqrt(Mito) ################### ##
## ######################################################################### ##

flog.info("Creating scatter plot...")

p6 <- ggplot(tmp, aes(x=umis, y=umis_mt_percent, color=genes)) +
  geom_point(size = 0.5) +
  scale_color_viridis() +
  theme_bw() +
  xlab("log10(Total UMIs)") +
  ylab("sqrt(Fraction of mitochondrial UMIs)") +
  # transform data
  scale_x_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
  scale_y_sqrt(breaks=c(0.00,0.05,0.20,0.25,0.50,0.75,1.00), labels=c(0.00,0.05,0.20,0.25,0.50,0.75,1.00)) +
  # add thresholds
  geom_hline(yintercept=c(0.2,0.5), linetype="dashed", color = "red", size=0.2) +
  geom_vline(xintercept=c(200,500,1000,5000), linetype="dashed", color = "red", size=0.2) +
  geom_hline(yintercept=c(opt$qc_minpercentmito, opt$qc_maxpercentmito), linetype="dashed", color = "blue", size=0.2) +
  # legend
  labs(color = "Number of genes") +
  theme(legend.position = "bottom",  legend.key.width = unit(1.5,"cm"), legend.key.height = unit(0.2,"cm")) +
  guides(color = guide_colorbar(title.position = "top"))
p6 <- ggMarginal(p6, type="density")

save_ggplots(file.path(opt$outdir, "log10umis_vs_sqrtMito"), p6)

#' log10(UMIs) vs sqrt(Fraction of mitochondrial UMIs) colored by Number Genes
#+ plot_log10umis_sqrtMito, include=TRUE, fig.height=5
print(p6)
#+ include=FALSE
#'

