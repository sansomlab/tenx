## Title ----
##
## Initial steps of the Seurat workflow for single-cell RNA-seq analysis
##
## Description ----
##
## This script performs the initial steps of the Seurat (http://satijalab.org/seurat/)
## workflow for single-cell RNA-seq analysis, including:
## (i) Reading in the data
## (ii) subsetting
## (iii) QC filtering
## (iv) Normalisation and removal of unwanted variation
## (v) Identification of variable genes
## (vi) Dimension reduction (PCA)
##
## The seurat object is saved as the R object "rds" in the output directory

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(xtable),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
        c("--tenxdir"),
        help="Location of the input 10x matrix"
        ),
    make_option(
        c("--project"),
        default="SeuratAnalysis",
        help="project name"
        ),
    make_option(
        c("--outdir"),
        default="seurat.out.dir",
        help="Location for outputs files. Must exist."
        ),
    make_option(
        c("--metadata"),
        default="none",
        help=paste(
            "A tab delimited text file containing cell metadata.",
            "A barcode column must be present and must match the barcodes of the 10x matrix"
            )
        ),
    # TODO
    # make_option(
    #     c("--batch"),
    #     default=F,
    #     help="If true, a single-column batch.tsv file must be present in the --tenxdir"
    #     ),
    make_option(
        c("--groupby"),
        default="none",
        help=paste(
            "The name of a column in the metadata table by which to group samples",
            "(e.g. sample_id)."
        )),
    make_option(
        c("--mingenes"),
        type="integer",
        default=200,
        help="min.genes"),
    make_option(
        c("--mincells"),
        type="integer",
        default=3,
        help=paste(
            "Include genes with detected expression in at least this many cells.",
            "See Seurat::CreateSeuratObject(min.cells=...)."
        )),
    make_option(
        c("--downsamplecells"),
        default=FALSE,
        help=paste(
            "Whether to randomly downsample cells by sample_id",
            "to match the sample with the lowest count of cells."
        )),
    #KRA: note that filtering _could_ be done before downsampling (with some obvious caveats)
    make_option(
        c("--qcmingenes"),
        type="integer",
        default=500,
        help=paste(
            "Minimal count of genes detected to retain a cell.",
            "See Seurat::FilterCells(subset.names='nGene', low.thresholds= ...)"
            )
        ),
    make_option(
        c("--qcmaxpercentmito"),
        type="double",
        default=0.05,
        help=paste(
            "Maximal percentage of UMI assigned to mitochondrial genes.",
            "See Seurat::FilterCells(subset.names='percent.mito', high.thresholds= ...)"
            )
        ),
    make_option(
        c("--latentvars"),
        default="nUMI,percent.mito",
        help=paste(
            "Latent variables to regress out.",
            "See Seurat::ScaleData(vars.to.regress=..., model.use=opt$modeluse)"
            )
        ),
    make_option(
        c("--modeluse"),
        default="linear",
        help=paste(
            "Model used to regress out latent variables.",
            "See Seurat::ScaleData(model.use=opt$modeluse)"
            )
        ),
    make_option(
        c("--sdcutoff"),
        type="double",
        default=0.5,
        help=paste(
            "Bottom cutoff on y-axis for identifying variable genes.",
            "See Seurat::FindVariableGenes(y.cutoff=...)"
            )
        ),
    make_option(
        c("--xlowcutoff"),
        type="double",
        default=0.1,
        help=paste(
            "Bottom cutoff on x-axis for identifying variable genes",
            "See Seurat::FindVariableGenes(x.low.cutoff=...)"
            )
        ),
    make_option(
        c("--xhighcutoff"),
        type="double",
        default=8,
        help=paste(
            "Top cutoff on x-axis for identifying variable genes",
            "See Seurat::FindVariableGenes(x.low.cutoff=...)"
            )
        ),
    make_option(
        c("--subsetcells"),
        default="use.all",
        help=paste(
            "A file containing the list of barcode ids to retain",
            "(no header, 1 per line)."
            )
        ),
    make_option(
        c("--subsetfactor"),
        default=NULL,
        help="A factor specified in metadata.tsv on which to subset"
        ),
    make_option(
        c("--subsetlevel"),
        default="none",
        help="The desired level of the sub-setting factor"),
    make_option(
        c("--cellcycle"),
        default="none",
        help="type of cell cycle regression to apply (none|all|difference)"
        ),
    make_option(
        c("--sgenes"),
        default="none",
        help=paste(
            "A vector of Ensembl gene ids associated with S phases.",
            "See Seurat::CellCycleScoring(s.genes=...)"
            )
        ),
    make_option(
        c("--g2mgenes"),
        default="none",
        help=paste(
            "A vector of Ensembl gene ids associated with G2M phase.",
            "See Seurat::CellCycleScoring(g2m.genes=...)"
            )
        )
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Input data ----

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##


#' Collect the number of cells through the workflow
#'
#' Initialises a data.frame of appends a new column
#' with a user-defined tag.
#'
#' @param s Seurat object
#' @param cell_numbers data.frame of an earlier call to
#' getCellNumbers, if applicable.
#' @param stage Character value to tag the new entry.
#'
#' @return A data.frame
getCellNumbers <- function(s, cell_numbers="none", stage="input") {
    ## cell_info <- getCellInfo(s)
    counts <- as.data.frame(table(s@meta.data$sample_id))
    colnames(counts) <- c("sample_id", stage)
    rownames(counts) <- counts$sample_id
    counts$sample_id <- NULL
    if ( identical(cell_numbers, "none") ) {
        result <- counts
        colnames(counts) <- stage
    } else {
        cnames <- c(colnames(cell_numbers), stage)
        result <- cbind(cell_numbers, counts[[stage]])
        colnames(result) <- cnames
    }
    return(result)
}

stats <- list()

cat("Importing matrix from: ", opt$tenxdir, " ... ")
data <- Read10X(opt$tenxdir)
cat("Done.\n")

## Seurat discards the Ensembl IDs and makes it's own identifiers (!)
## From https://github.com/satijalab/seurat/blob/master/R/preprocessing.R:
##
## rownames(data) <- make.unique(
##    names=as.character(
##        x=sapply(
##            X=gene.names,
##            FUN=ExtractField,
##            field=2,
##            delim="\\t"
##
## We need to track the seurat id -> Ensembl id mapping, e.g. for downstream GO analysis
## (and to guarentee reproducibility, e.g. between annotations versions)

inFile <- file.path(opt$tenxdir, "genes.tsv")
cat("Importing gene information from: ", inFile, " ... ")
genes <- read.table(inFile, as.is=TRUE)
cat("Done.\n")

colnames(genes) <- c("gene_id", "gene_name")
## TODO: consider scater::uniquifyFeatureNames()
genes$seurat_id <- as.character(make.unique(genes$gene_name))
rownames(genes) <- genes$seurat_id

write.table(
    genes, file.path(opt$outdir, "annotation.txt"),
    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
    )

## Initialize the Seurat object with the raw (non-normalized data)
## Note that this is slightly different than the older Seurat workflow,
## where log-normalized values were passed in directly.
## You can continue to pass in log-normalized values,
## just set do.logNormalize=F in the next step.
#s <- new("seurat", raw.data=data)

## Keep all genes expressed in >= 3 cells,
## keep all cells with >= 200 genes
## Perform log-normalization, first scaling each cell
## to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
cat("Creating Seurat object ... ")
s <- CreateSeuratObject(
    raw.data=data, min.cells=opt$mincells, min.genes=opt$mingenes,
    ## do.scale=F,
    ## do.center=F,
    #save.raw=T,
    project=opt$project
    )
cat("Done.\n")

# TODO: consider better alternatives (e.g. SummarizedExperiment has a rowData slot for that)
s@misc <- genes

## Read in the metadata
if ( identical(opt$metadata, "none") ) {
    stop("No metadata file given")
}

cat("Importing metadata ... ")
metadata <- read.table(opt$metadata, sep="\t", header=TRUE, as.is=TRUE)
cat("Done.\n")

rownames(metadata) <- metadata$barcode
metadata$barcode <- NULL

metadata <- metadata[s@cell.names, ]

# TODO: probably a more elegant way to show/handle this
if (!identical(rownames(metadata), s@cell.names)) {
    ## print out some debugging information
    cat("Count of rownames in metadata: ", length(rownames(metadata)), "\n")
    cat("Count of cell.names in Seurat object: ", length(s@cell.names), "\n")
    cat("--\n")
    cat("First rownames in metadata:\n")
    print(head(rownames(metadata)))
    cat("First cell.names in Seurat object:\n")
    print(head(s@cell.names))
    cat("--\n")
    cat("Last rownames in metadata:\n")
    print(tail(rownames(metadata)))
    cat("Last cell.names in Seurat object:\n")
    print(tail(s@cell.names))

    stop("Metadata barcode field does not match cell.names")
}

cat("Adding meta data ... ")
s <- AddMetaData(s, metadata)
cat("Done.\n")

## optionally add batch metadata
## if(opt$batch)
## {
##     batch <- read.table(file.path(opt$tenxdir, "batch.tsv"),
##                         header=FALSE, as.is=TRUE)$V1

##     barcodes <- read.table(file.path(opt$tenxdir, "barcodes.tsv"),
##                            header=FALSE, as.is=TRUE)$V1

##     names(batch) <- barcodes
##     print("Adding batch meta data")
##     print(head(batch))
##     s <- AddMetaData(s, batch, "batch")
## }

## ######################################################################### ##
## ####################### (ii) Subsetting ################################# ##
## ######################################################################### ##


# KRA: pipeline.yml says "to retain all the cells for a sample use the
# string "use.all".. why is it used for subsetting here?!
# TODO: use comments to explain what is done
if (opt$subsetcells!="use.all") {
    cells_to_retain <- scan(opt$subsetcells, "character")
    s <- SubsetData(s, cells.use=cells_to_retain)
} else {
    if (!is.null(opt$subsetfactor)) {
        if(!opt$subsetfactor %in% colnames(s@meta.data)) {
            stop("The given subsetting factor must match a column in the metadata")
        }

        if (!opt$subsetlevel %in% s@meta.data[[opt$subsetfactor]]) {
            stop("The specified level of the subsetting factor does not exist")
        }

        cells_to_retain <- rownames(s@meta.data)[
            s@meta.data[, opt$subsetfactor] == opt$subsetlevel
            ]

        if ( identical(length(cells_to_retain), 0L) ) {
            stop("No cells present in subset")
        }

        s <- SubsetData(s, cells.use=cells_to_retain)
    }
}

## ######################################################################### ##
## ##################### (iii) QC filtering ################################ ##
## ######################################################################### ##

## The number of genes and UMIs (nGene and nUMI) are
## automatically calculated for every object by Seurat.
## For non-UMI data, nUMI represents the sum of the
## non-normalized values within a cell
## We calculate the percentage of mitochondrial genes here
## and store it in percent.mito using the AddMetaData.
## The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
## NOTE: You must have the Matrix package loaded to calculate the percent.mito values.
## In humans the pattern is ^MT-, in mouse it is ^mt-: to accomodate both
## we simply ignore the case (specificity is still maintained).
mito.genes <- grep("^MT-", rownames(s@data), value=TRUE, ignore.case=TRUE)

## NOTE: unlike in the Seurat vignette, our data is not yet log transformed.
# TODO: is Matrix:: required here?
percent.mito <- Matrix::colSums(s@data[mito.genes, ]) / Matrix::colSums(s@data)

## AddMetaData adds columns to object@data.info,
## and is a great place to stash QC stats
s <- AddMetaData(s, percent.mito, "percent.mito")


if ( identical(opt$groupby, "none") ) {
    groupby <- NULL
} else {
    groupby <- opt$groupby
}

gp <- VlnPlot(
    s, c("nGene", "nUMI", "percent.mito"), size.title.use=14, size.x.use=12,
    group.by=groupby, x.lab.rot=TRUE, point.size.use=0.1, nCol=3
    )

# TODO: document where this function comes from
save_ggplots(
    file.path(opt$outdir, "qc.vlnPlot"), gp,
    width=7, height=4
    )

## GenePlot is typically used to visualize gene-gene relationships,
## but can be used for anything calculated by the object, i.e. columns
## in object@data.info, PC scores etc.
## Since there is a rare subset of cells with an outlier level of
## high mitochondrial percentage, and also low UMI content, we filter these as well

#' TODO: title
#'
#' TODO: describe
#'
#' @return TODO: describe
plot_fn <- function() {
    par(mfrow=c(1, 2))
    GenePlot(s, "nUMI", "percent.mito", cex.use=1)
    GenePlot(s, "nUMI", "nGene", cex.use=1)
    par(mfrow=c(1, 1))
}

# TODO: document where this function comes from
save_plots(
    file.path(opt$outdir, "qc.genePlot"), plot_fn=plot_fn,
    width=8, height=4
    )

# TODO: Document
cell_numbers <- getCellNumbers(s)

## We filter out cells that have unique gene counts over 2,500
## Note that accept.high and accept.low can be used to define a 'gate',
## and can filter cells not only based on nGene but on anything in the
## object (as in GenePlot above)

stats$no_cells <- ncol(s@data)
stats$qc_min_gene_threshold <- opt$qcmingenes
stats$qc_max_percent_mito_threshold <- opt$maxpercentmito

## s <- SubsetData(s, subset.name="nGene", accept.low=opt$qcmingenes)
## s <- SubsetData(s, subset.name="percent.mito", accept.high=opt$qcmaxpercentmito)

s <- FilterCells(
    s, subset.names=c("nGene", "percent.mito"),
    low.thresholds=c(opt$qcmingenes, -Inf),
    high.thresholds=c(Inf, opt$qcmaxpercentmito)
    )


stats$no_cells_after_qc <- ncol(s@data)

cell_numbers <- getCellNumbers(s, cell_numbers=cell_numbers, stage="after_qc_filters")

cat("Data dimensions after subsetting:\n")
print(dim(s@data))

## note that this overwrites pbmc@scale.data. Therefore, if you intend to use RegressOut,
## you can set do.scale=F and do.center=F in the original object to save some time.

# Optionally downsample cell numbers ----

if (as.logical(opt$downsamplecells)) {

    ##cell_info <- getCellInfo(s)

    if (!"sample_id" %in% colnames(s@meta.data)) {
        stop("sample_id must be specified to downsample cells")
        }

    mincells <- min(table(s@meta.data[[opt$groupby]]))

    cat(paste0("Downsampling to ", mincells, " per sample\n"))
    cells.to.use <- c()
    for (sample in unique(s@meta.data[[opt$groupby]])) {
        temp <- rownames(s@meta.data)[s@meta.data[[opt$groupby]] == sample]
        # KRA: We should probable set.seed() here, for reproducibility
        cells.to.use <- c(cells.to.use, sample(temp, mincells))
    }

    s <- SubsetData(s, cells.use=cells.to.use)

    ##cell_info <- getCellInfo(s)
    cell_numbers <- getCellNumbers(s, cell_numbers=cell_numbers, stage="after_downsampling")
    cat("Numbers of cells per sample after down-sampling:\n")
    print(cell_numbers)
}

print(
    xtable(cell_numbers, caption="Numbers of cells"),
    file=file.path(opt$outdir, "cell_numbers.tex")
    )

## ######################################################################### ##
## ######## (iv) Normalisation and removal of unwanted variation ########### ##
## ######################################################################### ##

s <- NormalizeData(
    object=s, normalization.method="LogNormalize", scale.factor=10E3
    )


latent.vars <- strsplit(opt$latentvars, ",")[[1]]
print(latent.vars)


# perform regression without cell cycle correction
s <- ScaleData(object=s, vars.to.regress=latent.vars,
               model.use=opt$modeluse)


## initialise the text snippet
tex = ""

## start building figure latex...
subsectionTitle <- getSubsubsectionTex("Cell cycle principal components analysis")
tex <- c(tex, subsectionTitle)


## If cell cycle genes are given, make a PCA of the cells based
## on expression of cell cycle genes

cell_cycle_genes <- FALSE

if (!is.null(opt$sgenes) & !is.null(opt$g2mgenes)){

  cell_cycle_genes <- TRUE
  cat("There are cell cycle genes")

  # get the genes representing the cell cycle phases
  sgenes_ensembl <- read.table(opt$sgenes, header=F, as.is=T)$V1
  sgenes <- s@misc$seurat_id[s@misc$gene_id %in% sgenes_ensembl]

  g2mgenes_ensembl <- read.table(opt$g2mgenes, header=F, as.is=T)$V1
  g2mgenes <- s@misc$seurat_id[s@misc$gene_id %in% g2mgenes_ensembl]

  # score the cell cycle phases
  s <- CellCycleScoring(
    object=s, s.genes=sgenes, g2m.genes=g2mgenes, set.ident=TRUE
  )

  s <- RunPCA(object = s, pc.genes = c(sgenes, g2mgenes), do.print = FALSE)

  ## PCA plot on cell cycle genes without regression
  gp <- PCAPlot(object = s)

  cc_plot_fn <- "cellcycle.without.regression.pca"
  cc_plot_path <- file.path(opt$outdir, cc_plot_fn)

  save_ggplots(cc_plot_path, gp, width=7, height=4)

  pcaCaption <- paste0("PCA analysis of cells based on expression of cell cycle genes ",
                       "prior to regression of cell-cyle effects")

  tex <- c(tex, getFigureTex(cc_plot_fn,
                             pcaCaption,
                             plot_dir_var="runDir"))

} else {
    tex <- c(tex, "Cell cycle genelists not supplied.")
}

## Perform regression (with or without cell cycle correction)
if ( identical(opt$cellcycle, "none") ) {

    cat("Data was scaled without correcting for cell cycle")

} else {

    if (!cell_cycle_genes){
        stop("Please provide lists of cell cycle sgenes and g2mgenes")
        }


    if ( identical(opt$cellcycle, "all") ){

      ## regress out all cell cycle effects
      s <- ScaleData(object=s, vars.to.regress=c(latent.vars, "S.Score", "G2M.Score"),
                     model.use=opt$modeluse)

        cat("Data was scaled to remove all cell cycle variation\n")

    } else if ( identical(opt$cellcycle, "difference") ) {

      ## regress out the difference between G2M and S phase scores
      s@meta.data$CC.Difference <- s@meta.data$S.Score - s@meta.data$G2M.Score
      s <- ScaleData(object=s, vars.to.regress=c(latent.vars, "CC.Difference"),
                     model.use=opt$modeluse)

      cat("Scaling was scaled to remove the difference between G2M and S phase scores\n")
    } else {
      stop("cellcycle regression type not understood")
    }

    ## visualise the cells by PCA of cell cycle genes after regression
    s <- RunPCA(object = s, pc.genes = c(sgenes, g2mgenes), do.print = FALSE)
    gp <- PCAPlot(object = s)

    cc_plot_fn <- paste("cellcycle.regressed", opt$cellcycle, "pca", sep=".")
    cc_plot_path <- file.path(opt$outdir, cc_plot_fn)

    save_ggplots(cc_plot_path, gp, width=7, height=4)

    pcaCaption <- paste0("PCA analysis of cells based on expression of cell cycle genes ",
                         "after regression of cell-cyle effects ",
                         "(regression type: ", opt$cellcycle, ")")

    tex <- c(tex, getFigureTex(cc_plot_fn,
                               pcaCaption,
                               plot_dir_var="runDir"))

}



tex_file <- file.path(opt$outdir, "cell.cycle.tex")

writeTex(tex_file, tex)



## ######################################################################### ##
## ################# (v) Identification of variable genes ################## ##
## ######################################################################### ##

## identify variable genes and output the plots
png(
    file.path(opt$outdir, "varGenesPlot.png"), width=8, height=5, unit="in",
    res=300
    )
s <- FindVariableGenes(
    s, mean.function=ExpMean, dispersion.function=LogVMR,
    x.low.cutoff=opt$xlowcutoff, x.high.cutoff=opt$xhighcutoff,
    y.cutoff=opt$sdcutoff, plot.both=TRUE, do.text=FALSE, do.contour=FALSE
    )
dev.off()

cat("no. variable genes: ", length(s@var.genes), "\n")
stats$no_variable_genes <- length(s@var.genes)

print(stats)

# write out some statistics into a latex table.
print(
    xtable(t(data.frame(stats)), caption="Run statistics"),
    file=file.path(opt$outdir, "stats.tex")
    )

## ######################################################################### ##
## ################### (vi) Dimension reduction (PCA) ###################### ##
## ######################################################################### ##

# perform PCA using the variable genes
s <- RunPCA(s, pc.genes=s@var.genes, pcs.compute=30, do.print=FALSE)

n_cells_pca <- min(1000, length(s@cell.names))

# Write out heatmaps showing the genes loading the first 12 components
plot_fn <- function() {
    PCHeatmap(
        s, pc.use=1:12, cells.use=n_cells_pca, do.balanced=TRUE,
        label.columns=FALSE, cexRow=0.8, use.full=FALSE
        )
}

save_plots(
    file.path(opt$outdir, "pcaComponents"), plot_fn=plot_fn,
    width=8, height=12
    )


# Write out the PCA elbow (scree) plot
png(
    file.path(opt$outdir, "pcaElbow.png"), width=5, height=4, units="in",
    res=300
    )
PCElbowPlot(s)
dev.off()

# Save the R object
saveRDS(s, file=file.path(opt$outdir, "begin.rds"))

message("Completed")
