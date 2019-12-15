#' ---
#' title: "Subsetting, QC and filtering"
#' output: 
#'  html_document:
#'   self_contained: false
#' params:
#'  task_yml: "/gfs/devel/ssansom/tenx/Rmd/seurat_quality_assessment.test.yml"
#'  fig_path: "fig.dir/"
#' ---
#' ---
#' Perform quality assessment of a single-cell RNA-seq dataset
#'
#' This script can be compiled from the command line or run interactively in
#' Rstudio.
#' ---


#+ setup, include=FALSE, echo=FALSE

# Libraries --------------------------------------------------------------------
stopifnot(require(yaml))
stopifnot(require(Seurat))
stopifnot(require(sctransform))
stopifnot(require(ggplot2))
stopifnot(require(dplyr))
stopifnot(require(Matrix))
stopifnot(require(tenxutils))
stopifnot(require(reshape2))
stopifnot(require(future))
stopifnot(require(gridExtra))
stopifnot(require(knitr))

# set chunk options
opts_chunk$set(echo=FALSE,
               warning=FALSE,
               message = FALSE,
               fig.path = params$fig_path)

# Parameters -------------------------------------------------------------------
# The script expects the following paramters:

default_options <- list(
  #project
  "project" = "single cell project",

  # Location of the input 10x matrix
  "tenxdir" = "",

  # location for the ouput files
  "outdir" = "",

  # Path to file containing the cell metadata ("barcode" column must be present)
  "metadata" = NULL,

  # The name of a column in the metadata table by which to group samples
  # (e.g. sample_id)                   #
  "qc_groupby" = "",

  "qc_initial_mingenes" = 200,

  # Include genes with detected expression in at least this many cells.
  # See Seurat::CreateSeuratObject(min.cells=...)."
  "qc_mincells" = 3,

  # Whether to randomly downsample cells so that each group has the
  # same number of cells."
  "downsamplecells" = FALSE,

  # The seed (an integer) to use when down sampling cells
  "qc_seed" = 42,

  # Minimal count of genes detected to retain a cell.
  # See Seurat::FilterCells(subset.names='nGene', low.thresholds= ...)
  "qc_mingenes" = 500,

  # Minimum percentage of UMI assigned to mitochondrial genes.
  # See Seurat::FilterCells(subset.names='percent.mito', low.thresholds= ...)
  "qc_minpercentmito" = -Inf,

  # Maximal percentage of UMI assigned to mitochondrial genes.",
  # See Seurat::FilterCells(subset.names = 'percent.mito', high.thresholds= ...)
  "qc_maxpercentmito" = 0.05,

  # Max nCount_RNA to retain a cell
  "qcmaxcount" = NULL,

  # The normalization method to use
  "normalization_method" = "log-normalization",

  # Latent variables to regress out.,
  # See Seurat::ScaleData(vars.to.regress=..., model.use=opt$modeluse)
  "regress_latentvars" = "nUMI,percent.mito",

  # Model used to regress out latent variables.,
  # See Seurat::ScaleData(model.use=opt$regress_modeluse)
  "regress_modeluse" = "linear",

  # Method for variable gene selection.
  # Either "top.genes" or "mean.var.plot"
  "vargenes_method" = "mean.var.plot",

  # Number of highly variable genes to retain
  "vargenes_topgenes" = 1000,

  # Bottom cutoff on y-axis for identifying variable genes.
  # See Seurat::FindVariableFeatures(y.cutoff=...)
  "vargenes_sdcutoff" = 0.5,

  # Bottom cutoff on x-axis for identifying variable genes
  # See Seurat::FindVariableFeatures(x.low.cutoff=...)
  "vargenes_xlowcutoff" = 0.1,

  # Top cutoff on x-axis for identifying variable genes
  # See Seurat::FindVariableFeatures(x.low.cutoff=...)
  "vargenes_xhighcutoff" = 8,

  # minimum mean of log counts when using trendvar method
  "vargenes_minmean" = 0,

  # significance threshold for trendvar method
  "vargenes_padjust" = 0.05,

  # Subset to a whitelist.
  # Specify the path to a file containing a list of barcode ids
  "subset_white_list" = NULL,

  # Subset to a blacklist
  # Specify the path to a file containing a list of barcode ids
  "subset_black_list" = NULL,

  # A factor specified in metadata.tsv on which to subset
  "subset_factor" = NULL,

  # The desired level of the sub-setting factor
  "subset_factor_level" = "none",

  # A vector of Ensembl gene ids associated with S phases.
  # See Seurat::CellCycleScoring(s.genes=...)
  "cellcycle_sgenes" = "none",

  # A vector of Ensembl gene ids associated with G2M phase.
  # See Seurat::CellCycleScoring(g2m.genes=...)
  "cellcycle_g2mgenes" = "none",

  # Number of cores to be used for the Jackstraw analysis
  "numcores" = 12

)

# If running interactively comment this line out or specify the
# location of a file containing the parameters that you wish to use, e.g.
# options <- read_yaml("/gfs/devel/ssansom/tenx/Rmd/seurat_quality_assessment.test.yml")

options <- read_yaml(params$task_yml)

# if(! all(names(options) %in% names(default_options))) {
#  stop("option not recognised")
# }

# Update the default options
if(!is.null(options)) {
  opt <- utils::modifyList(default_options, options)
} else{
  opt <- default_options
}

message("Running with options:\n")
print(opt)

# Set up for multiprocessing
plan("multiprocess",
     workers = opt$numcores)


# Function definitions ---------------------------------------------------------

.qc_plots <- function(s, group_by = "") {
  # Draw some violin plots
  gp <- VlnPlot(s,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
          group.by=group_by,
          pt.size = 0.1,
          ncol=3)

  print(gp)

  # draw some scatter plots
  gp1 <- FeatureScatter(s, "nCount_RNA", "percent.mito", pt.size=1)
  gp2 <- FeatureScatter(s, "nCount_RNA", "nFeature_RNA", pt.size=1)
  grid.arrange(grobs=list(gp1+ theme(legend.position="none"),
                          gp2+ theme(legend.position="none")), ncol=2)
}

is.not.set <- function(x="none")
{
    if(is.null(x)){
      return(TRUE)
    } else if (x == "none" | x == "NULL") {
      return(TRUE)
    } else {
      return(FALSE)
    }
}

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##

stats <- list()

data <- Read10X(opt$tenxdir)

s <- CreateSeuratObject(counts=data,
                        min.cells=opt$qc_mincells,
                        min.features=opt$qc_initial_mingenes,
                        project=opt$project)

# Add gene information
features <- read.table(file.path(opt$tenxdir, "features.tsv.gz"),
                    as.is=TRUE)

colnames(features) <- c("gene_id", "gene_name")
features$seurat_id <- as.character(make.unique(features$gene_name))
rownames(features) <- features$seurat_id
s@misc <- features

## Read in the metadata
if ( is.not.set(opt$metadata) ) {
    stop("No metadata file given")
}

s <- seurat_add_metadata(s, opt$metadata)

# ensure that the grouping factor is present in the metadata
if (!opt$qc_groupby %in% colnames(s[[]])) {
    stop(paste("The specified grouping factor:",
               opt$qc_groupby,
               "was not found in the metadata",
               sep=" "))
}

cell_numbers <- seurat_track_cell_numbers(s, groupby=opt$qc_groupby)

## ######################################################################### ##
## ####################### (ii) Subsetting ################################# ##
## ######################################################################### ##

#' ## Data subsetting
#'
#' The options for data subsetting were:
#'
#' - whitelist: `r opt$subset_white_list`
#' - blacklist: `r opt$subset_black_list`
#' - subset factor: `r opt$subset_black_list`
#'
#' Before subsetting there were: `r length(Cells(s))` cells.

if(!is.not.set(opt$subset_white_list)
  | !is.not.set(opt$subset_black_list)
  | !is.not.set(opt$subset_factor)){
  s <- seurat_subset_cells(s,
                         subset_white_list = opt$white_list,
                         subset_black_list = opt$black_list,
                         subset_factor = opt$subset_factor,
                         subset_factor_level = opt$subset_factor_level
                         )
}

#' After subsetting theere were `r length(Cells(s))` cells.
cell_numbers <- seurat_track_cell_numbers(s,
                                          cell_numbers=cell_numbers,
                                          stage="after_subsetting",
                                          groupby=opt$qc_groupby)

## ######################################################################### ##
## ################# (iii) Pre filtering QC plots ########################## ##
## ######################################################################### ##

#' ## Pre-filtering quality assessment
#'
s <- seurat_pc_mito(s)

#+ fig.height=3
.qc_plots(s, group_by=opt$qc_groupby)



## ######################################################################### ##
## ######### (iv) Data Filtering and downsampling of cell numbers ########## ##
## ######################################################################### ##

#' ## Filtering and downsampling
#'
#' The specified filtering paramters were:
#'
#' - minimum no. genes: `r opt$qc_mingenes`
#' - mininimum % mitochondrial genes: `r opt$qc_minpercentmito`
#' - maximum % mitochondrial genes: `r opt$qc_maxpercentmito`
#' - maximum no. RNA molecules: `r opt$qcmaxcount`
#'
#' Before filtering there were: `r length(Cells(s))` cells.

if (! is.not.set(opt$qcmaxcount)) {
  s <- subset(s, subset = nFeature_RNA > opt$qc_mingenes &
                   percent.mito > opt$qc_minpercentmito &
                   percent.mito < opt$qc_maxpercentmito &
                   nCount_RNA < opt$qcmaxcount)
  stats$qc_max_count_threshold <- opt$qcmaxcount
  } else {
    s <- subset(s, subset = nFeature_RNA > opt$qc_mingenes &
                     percent.mito > opt$qc_minpercentmito &
                     percent.mito < opt$qc_maxpercentmito)
    }

#' After filtering there were: `r length(Cells(s))` cells.


stats$no_cells_after_qc <- ncol(GetAssayData(object = s))

cell_numbers <- seurat_track_cell_numbers(s,
                               cell_numbers=cell_numbers,
                               stage="after_qc",
	     		       groupby=opt$qc_groupby)

# Optionally downsample cell numbers ----

if(is.not.set(opt$qc_seed))
{
    qc_seed <- sample(1:2^15,1)
    set.seed(qc_seed)
    message("Seed set to: ", qc_seed)
} else {
    qc_seed <- opt$qc_seed
}

if (as.logical(opt$downsamplecells)) {

    mincells <- min(table(s[[]][,opt$qc_groupby]))

    cat(paste0("Downsampling to ", mincells, " per sample\n"))
    cells.to.use <- c()
    for (group in unique(s[[]][,opt$qc_groupby])) {
        print(group)
        temp <- rownames(s[[]])[s[[]][,opt$qc_groupby] == group]
        print(head(temp))

        cells.to.use <- c(cells.to.use, sample(temp, mincells))
        print(length(cells.to.use))
    }

    s <- SubsetData(s,
                    cells=cells.to.use)

    cell_numbers <- getCellNumbers(s, cell_numbers=cell_numbers,
    		    		   stage="after_downsampling", groupby=opt$qc_groupby)

    cat("Numbers of cells per sample after down-sampling:\n")
    print(cell_numbers)

    cell_numbers <- seurat_track_cell_numbers(s,
                                              cell_numbers=cell_numbers,
                                              stage="after_downsampling",
                                              groupby=opt$qc_groupby)
}


# Save a list of the selected cells
# (after subsetting, qc and downsampling)
fileCon <- gzfile(file.path(opt$outdir, "selected.cells.txt.gz"))
write.table(Cells(s),fileCon)

## ######################################################################### ##
## ######################## (v) Post filtering QC plots #################### ##
## ######################################################################### ##

#' ## Post filtering (and down-sampling) quality assessment
#'
#+ fig.height=3
.qc_plots(s, group_by = opt$qc_groupby)

#' ## Summary of cell numbers
#'
kable(cell_numbers, caption="Numbers of cells")

## ######################################################################### ##
## ###### (vi) Visualise cell cycle effects in the cells passing QC ######### ##
## ######################################################################### ##

#' ## Visualisation of cell cycle effects
#'
#' The following lists of cell cycle genes were given:
#'
#' - S phase genes: `r opt$cellcycle_sgenes`
#' - G2M phase genes: `r opt$cellcycle_g2mgenes`
#'
#' The plot below visualises cell cycle effects before regression for
#' cell cycle genes.

cell_cycle_genes <- FALSE

#+ fig.height=3, fig.width=4
if (!(is.not.set(opt$cellcycle_sgenes) | opt$cellcycle_sgenes=="none")
    & !(is.not.set(opt$cellcycle_g2mgenes) | opt$cellcycle_g2mgenes=="none"))
{
  latent.vars <- strsplit(opt$regress_latentvars, ",")[[1]]
  print(latent.vars)

  ## Currently, we need to log-normlize and scale the RNA assay
  ## for gene-level analyses even if we use sctransform
  ## for characterisation of the cells/subpopulations.
  ##
  ## see e.g. https://github.com/satijalab/seurat/issues/1717
  ## and  https://github.com/satijalab/seurat/issues/1421


  message("Performing initial log-normalization")

  ## Perform log-normalization of the RNA assay
  s <- NormalizeData(object=s,
                     normalization.method="LogNormalize",
                     scale.factor=10E3)

  ## Initial scaling of the RNA assay data
  all.genes <- rownames(s)
  s <- ScaleData(object=s,
                 features = all.genes,
                 vars.to.regress=latent.vars,
                 model.use=opt$regress_modeluse)

  VariableFeatures(object = s) <- rownames(s)


  cell_cycle_genes <- TRUE

  # get the genes representing the cell cycle phases
  sgenes_ensembl <- read.table(opt$cellcycle_sgenes, header=F, as.is=T)$V1
  sgenes <- s@misc$seurat_id[s@misc$gene_id %in% sgenes_ensembl]

  g2mgenes_ensembl <- read.table(opt$cellcycle_g2mgenes, header=F, as.is=T)$V1
  g2mgenes <- s@misc$seurat_id[s@misc$gene_id %in% g2mgenes_ensembl]

  # score the cell cycle phases
    s <- CellCycleScoring(object=s,
                          s.features=sgenes,
                          g2m.features=g2mgenes,
                          set.ident=TRUE)

    s <- RunPCA(object = s,
                features = c(sgenes, g2mgenes),
                do.print = FALSE)

  ## PCA plot on cell cycle genes without regression
  gp <- PCAPlot(object = s)
  gp
}

message("Completed")
