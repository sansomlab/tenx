## Characterise a set of differentially expressed genes (per cluster)
## This script will produce:
## (1) MA plot (colored by p-value)
## (2) MA style plot for frequency of expression (colored by p-value)
## (3) Fold change vs Frequency change plot (coloured by p-value)
## (4) Traditional volcano plot
## (5) Violin plots of the top DE genes.
##
## The idea is to allow evaluation of the performance of the
## DE alogorithm applied (to aid/inform alogorithm selection)
## run only specified comparison (in order for parallel execution)

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(gplots),
  require(reshape2),
  require(xtable),
  require(ggplot2),
  require(ggrepel),
  require(gridExtra),
  require(dplyr),
  require(tenxutils)
)


# Options ----

option_list <- list(
    make_option(c("--degenes"), default="begin.Robj",
                help="Summary table of differentially expressed genes"),
    make_option(c("--seuratobject"), default="begin.Robj",
                help="seurat object"),
    make_option(c("--seuratassay"), default="RNA",
                help="the assay to set as default (used for the violin plots)"),
    make_option(c("--clusterids"), default="begin.Robj",
                help="clusterids"),
    make_option(c("--cluster"), default="none",
                help="The cluster to characterise"),
    make_option(c("--testfactor"), default="NULL",
                help="a metadata factor used to group the violin plots"),
    make_option(c("--a"), default="NULL",
                help="first level of constrast"),
    make_option(c("--b"), default="NULL",
                help="second level of constrast"),
    make_option(c("--plotdirvar"), default="clusterMarkerDEPlotsDir",
                help="latex var containing location of the plots"),
    make_option(c("--useminfc"), default=FALSE,
                help="Use minimum fold change for second page of violin plots"),
    make_option(c("--ncol"), type="integer", default=4,
                help="Number of columns for the plots"),
        make_option(c("--nrow"), type="integer", default=4,
                help="Number of rows for the plots"),
    make_option(c("--pointsize"), type="integer", default=FALSE,
                help="pointsize for the violin plots"),
    make_option(c("--pdf"), default=FALSE,
                help="Produce pdf plots"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$testfactor=="NULL") { opt$testfactor <- NULL }

cat("Running with options:\n")
print(opt)

# set the run specs
run_specs <- paste(opt$numpcs,opt$resolution,opt$algorithm,opt$testuse,sep="_")

## set up the ggplot theme
theme_set(theme_classic(base_size = 10))

if(is.null(opt$testfactor))
{
    file_suffix = NULL
    a = "cluster_mean"
    b = "other_mean"
} else {
    file_suffix = paste0("between")
    a = paste(opt$a,"mean",sep="_")
    b = paste(opt$b,"mean",sep="_")
}

## initialise the text snippet
tex = ""

## read in the de genes

cluster <- opt$cluster
data <- read.table(gzfile(opt$degenes),header=T,as.is=T,sep="\t")

if(!cluster %in% data$cluster)
{
    ## we have no DE genes.
    message("No significant genes for cluster: ",cluster)
} else

{

data <- data[data$cluster == cluster,]


message("processing cluster: ",cluster)

message("making scatter plots")
## analyse expression level vs fold change
## note that we intentionally ignore the pseudocount that seurat has added.

#data <- decluster
print(head(data))

data$M <- log2( exp(data[[a]]) / exp(data[[b]]) )
data$A <- 1/2 * (log2(exp(data[[a]])) + log2(exp(data[[b]])))

ma <- plotMA(data,xlab="expression level (log2)",ylab="fold change (log2)")
vol <- plotVolcano(data,xlab="fold change (log2)",ylab="adjusted p-value (-log10)")

data$M <- NULL
data$A <- NULL

## analyse percentage vs change in percentage
#data <- decluster

data$M <- log2((data$pct.1*100 + 1) / (data$pct.2*100 + 1))
data$A <- 1/2 * (log2(data$pct.1*100 + 1) + log2(data$pct.2*100 + 1))

pma <- plotMA(data,xlab="percent cells (log2)",ylab="percent change (log2)")
pvol <- plotVolcano(data,xlab="percent change (log2)",ylab="adjusted p-value (-log10)")

data$M <- NULL
data$A <- NULL

## directly compare fold change and percentage change
#data <- decluster
data$L <- log2((data$pct.1*100 + 1) / (data$pct.2*100 + 1))
data$M <- log2(exp(data[[a]]) / exp(data[[b]]))

fve <- plotFvE(data)

data$M <- NULL
data$L <- NULL

## sneakily extract the legend...
legend <- g_legend(fve)
fve <- fve + theme(legend.position = 'none')

gs <- list(ma, vol, pma, pvol, fve, legend)
rm(ma, vol, pma, pvol, fve, legend)

lay <- rbind(c(1,2),c(3,4),c(5,6))
ga <- arrangeGrob(grobs = gs, layout_matrix = lay)
rm(gs)

defn <- paste(c("dePlots",cluster,file_suffix),collapse=".")
defpath <- file.path(opt$outdir, defn)

save_ggplots(defpath,
             ga,
             to_pdf=opt$pdf,
             width=8,
             height=10)

rm(ga)
gc()

## start building figure latex...
subsectionTitle <- getSubsectionTex(paste0("Cluster ",cluster,": summary plots"))
tex <- c(tex, subsectionTitle)

deCaption <- paste("Differential expression summary plots for cluster ",cluster)
tex <- c(tex, getFigureTex(defn, deCaption,plot_dir_var=opt$plotdirvar))

## make two pages of violin plots for the top average cluster markers:
## Page (1) (a) top +ve by p-value
##          (b) top +ve by fold change (non redundant with (a))
## Page (2) (a) top -ve by p-value
##          (b) top -ve by fold change (non redundant with (b))
##
## - Only significant genes will be plotted
## - For cluster marker genes, min fold change (vs all other clusters is used).

## order by should be one of "p-value", "fold-change" or "min-fold-change".
## enforce significance

## read in the seurat object
s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
Idents(s) <- cluster_ids

## set the default assay
message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("seurat_characteriseClusterDEGenes.R running with default assay: ", DefaultAssay(s))

violin_fn_template <- paste(c("violinPlots.type",cluster,file_suffix),collapse=".")

data <- data[data$p.adj < 0.05,]

print(dim(data))

if(opt$useminfc)
{
    fc_type="minimum fold-change vs all other clusters"
} else {
    fc_type="fold change"
}

ncol = opt$ncol

if(is.null(opt$testfactor))
{
    analysis_title = " marker genes"
    ident.include <- NULL
} else {
    analysis_title = "ly differentially expressed genes"
    ident.include <- cluster
}

## make the +ve plots
message("making violin plots for +ve genes")
pos_tex <- violinPlotSection(data, s, cluster_ids, type="positive",
                             group.by = opt$testfactor,
                             ident.include = ident.include, vncol=opt$ncol, vnrow=opt$nrow,
                             pt_size=opt$pointsize,
                             outdir = opt$outdir,
                             analysis_title = analysis_title, fc_type = fc_type,
                             use.minfc = opt$useminfc,
                             plot_dir_var=opt$plotdirvar,
                             to_pdf=opt$pdf)

## make the -ve plots
message("making violin plots for -ve genes")
neg_tex <- violinPlotSection(data, s, cluster_ids, type="negative",
                             group.by = opt$testfactor,
                             ident.include = ident.include, vncol=opt$ncol, vnrow=opt$nrow,
                             pt_size=opt$pointsize,
                             outdir = opt$outdir,
                             analysis_title = analysis_title, fc_type = fc_type,
                             use.minfc = opt$useminfc,
                             plot_dir_var=opt$plotdirvar,
                             to_pdf=opt$pdf)

tex <- c(tex, pos_tex, neg_tex)

tex_file <- file.path(opt$outdir,
                      paste(c("characterise.degenes",cluster,file_suffix,"tex"),collapse="."))

writeTex(tex_file, tex)

}

message("seurat_characteriseClusterDEGenes.R final default assay: ", DefaultAssay(s))

message("completed")
