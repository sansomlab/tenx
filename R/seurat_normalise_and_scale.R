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
  require(sctransform),
  require(ggplot2),
  require(dplyr),
  require(Matrix),
  require(xtable),
  require(tenxutils),
  require(reshape2),
  require(future)
)


# Options ----

option_list <- list(
    make_option(
        c("--seuratobject"),
        help="Location of the seuratobject"
    ),
    make_option(
        c("--project"),
        default="SeuratAnalysis",
        help="project name"
    ),
    make_option(
        c("--species"),
        default="not_set",
        help="species, e.g. mm or hs"
    ),
    make_option(
        c("--outdir"),
        default="seurat.out.dir",
        help="Location for outputs files. Must exist."
        ),
    make_option(
        c("--seed"),
        default=NULL,
        help=paste(
            "The seed (an integer) to use when down sampling the cells"
        )),
    make_option(
        c("--normalizationmethod"),
        default="log-normalization",
        help="The normlization method to use"
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
        c("--vargenesmethod"),
        default="vst",
        help=paste(
            "Method for variable gene selection.",
            "Either top.genes or mean.var.plot"
            )
    ),
    make_option(
            c("--regressgenes"),
        default="HVG",
        help=paste(
            "Genes to use for scaling the data either all or HVG"
            )
        ),
    make_option(
        c("--topgenes"),
        type="integer",
        default=1000,
        help=paste(
            "Number of highly variable genes to retain"
            )
        ),
    make_option(
        c("--sdcutoff"),
        type="double",
        default=0.5,
        help=paste(
            "Bottom cutoff on y-axis for identifying variable genes.",
            "See Seurat::FindVariableFeatures(y.cutoff=...)"
            )
        ),
    make_option(
        c("--xlowcutoff"),
        type="double",
        default=0.1,
        help=paste(
            "Bottom cutoff on x-axis for identifying variable genes",
            "See Seurat::FindVariableFeatures(x.low.cutoff=...)"
            )
        ),
    make_option(
        c("--xhighcutoff"),
        type="double",
        default=8,
        help=paste(
            "Top cutoff on x-axis for identifying variable genes",
            "See Seurat::FindVariableFeatures(x.low.cutoff=...)"
            )
    ),
        make_option(
        c("--minmean"),
        type="double",
        default=0,
        help=paste(
            "minimum mean of log counts when using trendvar method"
            )
        ),
        make_option(
        c("--vargenespadjust"),
        type="double",
        default=0.05,
        help=paste(
            "significance threshold for trendvar method"
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
        c("--cellcycle"),
        default="none",
        help="type of cell cycle regression to apply (none|all|difference)"
        ),
    make_option(
        c("--sgenes"),
        default=NULL,
        help=paste(
            "A vector of Ensembl gene ids associated with S phases.",
            "See Seurat::CellCycleScoring(s.genes=...)"
            )
        ),
    make_option(
        c("--g2mgenes"),
        default=NULL,
        help=paste(
            "A vector of Ensembl gene ids associated with G2M phase.",
            "See Seurat::CellCycleScoring(g2m.genes=...)"
            )
    ),
    make_option(
        c("--numcores"),
        type="integer",
        default=12,
        help="Number of cores to be used for the Jackstraw analysis"
    ),
    make_option(
        c("--memory"),
        type="integer",
        default=4000,
        help="Amount of memory (mb) to request"
    ),

    make_option(c("--usesigcomponents"), default=FALSE,
                help="use significant principle component"),
    make_option(c("--components"), type="integer", default=10,
                help="if usesigcomponents is FALSE, the number of principle components to use"),
    make_option(
        c("--plotdirvar"),
        default="sampleDir",
        help="latex var containig plot location"
    )


    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

if(!opt$normalizationmethod %in% c("log-normalization", "sctransform"))
{
    stop(paste("normalisation method not recognised. Supported methods are",
               "`log-normalization` and `sctransform`"))
}

plan("multiprocess",
     workers = opt$numcores)

#plan("sequential")

options(future.globals.maxSize = opt$memory * 1024^2)

# Input data ----

s <- readRDS(opt$seuratobject)

## ######################################################################### ##
## # (i) Initial normalisation, variable gene identification and scaling ## ##
## ######################################################################### ##


## No cell cycle correction is applied at this stage.

if(opt$latentvars == "none") {
    latent.vars = NULL
    message("no latent vars specified")
} else {
    latent.vars <- strsplit(opt$latentvars, ",")[[1]]
    message("latent vars: ", latent.vars)
}

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

message("log-normalisation completed")

gc()


if (opt$cellcycle != "none")
{
    ## see https://github.com/satijalab/seurat/issues/1679
    ## requires normalized data and should be done prior to
    ## sctransform
    cc_genes <- getCellCycleGenes(sgenes_file = opt$sgenes,
                              g2mgenes_file = opt$g2mgenes,
                              species = opt$species)

    ## score the cell cycle phases
    s <- CellCycleScoring(object=s,
                          s.features=cc_genes$s.genes,
                          g2m.features=cc_genes$g2m.genes,
                          set.ident=TRUE)
 }


if (opt$cellcycle == "none")
{
    vars.to.regress <- c(latent.vars)
} else if ( identical(opt$cellcycle, "all") )
{

    message("Cell cycle correction for S and G2M scores will be applied\n")
    vars.to.regress <- c(latent.vars, "S.Score", "G2M.Score")

} else if ( identical(opt$cellcycle, "difference") ) {

    message("Cell cycle correction for the difference between G2M and S phase scores will be applied\n")
    s$CC.Difference <- s$S.Score - s$G2M.Score
    vars.to.regress <- c(latent.vars, "CC.Difference")

} else {
    stop("Cell cycle regression type not recognised")
}


## ########################################################################### #
## ###################### Normalisation and scaling ########################## #
## ########################################################################### #

## log-normalisation and scaling is always performed as we use the RNA slot
## for differential expression and visualisation downstream.

s <- FindVariableFeaturesMod(seurat_object=s,
                             method=opt$vargenesmethod,
                             nfeatures=opt$topgenes,
                             xlowcutoff=opt$xlowcutoff,
                             xhighcutoff=opt$xhighcutoff,
                             sdcutoff=opt$sdcutoff,
                             minmean=opt$minmean,
                             padj=opt$vargenespadjust)


message("scaling the data")
    ## Initial scaling of the RNA assay data

## variable gene identification
if(opt$regressgenes=="all") {
    feature.genes <- rownames(s)
} else {
    feature.genes <- NULL }

message("Performing log-normalization")

s <- ScaleData(object=s,
               features = feature.genes,  # if NULL defaults to HVG
               vars.to.regress=vars.to.regress,
               model.use=opt$modeluse,
               assay="RNA")


if(opt$normalizationmethod=="sctransform")
{
    message("Performing SCTransform normalization")

    s <- SCTransform(object=s,
                     assay="RNA",
                     new.assay.name="SCT",
                     do.correct.umi=TRUE,
                     variable.features.n=3000,
                     vars.to.regress=vars.to.regress,
                     do.scale=FALSE,
                     do.center=TRUE,
                     return.only.var.genes=FALSE)

    ## Note that the SCT slot will now be set as default.

}


message("seurat_begin.R object final default assay: ", DefaultAssay(s))

# Save the R object
saveRDS(s, file=file.path(opt$outdir, "begin.rds"))

message("Completed")
