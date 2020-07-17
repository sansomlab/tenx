## Title ----
##
## This script makes plots of the HVG and,  optionally, of variation due to cell
## cycle genes
##
## Description ----
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
        c("--species"),
        default="not_set",
        help="species, e.g. mm or hs"
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
        c("--cellcycle"),
        default="none",
        help="type of cell cycle regression to apply (none|all|difference)"
        ),
    make_option(
        c("--sgenes"),
        default= NULL,
        help=paste(
            "A vector of Ensembl gene ids associated with S phases.",
            "See Seurat::CellCycleScoring(s.genes=...)"
            )
        ),
    make_option(
        c("--g2mgenes"),
        default= NULL,
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

plan("multiprocess",
     workers = opt$numcores)

#plan("sequential")

options(future.globals.maxSize = opt$memory * 1024^2)

# Input data ----

# Input data ----

if (endsWith(opt$seuratobject, ".rds")) {
  message(sprintf("readRDS: %s", opt$seuratobject))
  s <- readRDS(opt$seuratobject)
} else {
  message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
  stopifnot(require(SeuratDisk))
  s <- LoadH5Seurat(opt$seuratobject)
}


## ######################################################################### ##
## # (i) Initial normalisation, variable gene identification and scaling ### ##
## ######################################################################### ##

## No cell cycle correction is applied at this stage.
if(opt$latentvars == "none") {
    latent.vars = NULL
    message("no latent vars specified")
} else {
    latent.vars <- strsplit(opt$latentvars, ",")[[1]]
    message("latent vars: ", latent.vars)
}


if(opt$normalizationmethod=="log-normalization")
{

    message("Performing initial log-normalization")
    ## Perform log-normalization of the RNA assay
    s <- NormalizeData(object=s,
                   normalization.method="LogNormalize",
                   scale.factor=10E3)

    message("log-normalisation completed")

    s <- FindVariableFeaturesMod(seurat_object=s,
                                 method=opt$vargenesmethod,
                                 nfeatures=opt$topgenes,
                                 xlowcutoff=opt$xlowcutoff,
                                 xhighcutoff=opt$xhighcutoff,
                                 sdcutoff=opt$sdcutoff,
                                 minmean=opt$minmean,
                                 padj=opt$vargenespadjust)


    xthreshold <- opt$xlowcutoff
    if(opt$vargenesmethod=="trendvar")
    {
        xthreshold <- opt$minmean
    }

    message("scaling the data")

    ## Initial scaling of the RNA assay data
    ## For this purpose (cell cycle investigation we always use the HVG
    s <- ScaleData(object=s,
                   vars.to.regress=latent.vars,
                   model.use=opt$modeluse)

    ## Normalisation specific options.
    message("scaling complete")

} else if(opt$normalizationmethod=="sctransform")
{
    message("Performing initial SCTransform normalization")

    s <- SCTransform(object=s,
                     assay="RNA",
                     new.assay.name="SCT",
                     do.correct.umi=TRUE,
                     variable.features.n=3000,
                     vars.to.regress=latent.vars,
                     do.scale=FALSE,
                     do.center=TRUE,
                     return.only.var.genes=FALSE)

    ## Note that the SCT slot will now be set as default.

} else {
    stop("Invalid normalization method specified")
}


## make a plot that shows the variable genes.

xx <- HVFInfo(object = s)

xx$var.gene = FALSE
xx$var.gene[rownames(xx) %in% VariableFeatures(object = s)] <- TRUE

if(opt$normalizationmethod=="log-normalization")
{
    if(opt$vargenesmethod=="vst") {
        yvar = "variance.standardized"
    } else { yvar = "dispersion" }

    xvar = "mean"

} else if(opt$normalizationmethod=="sctransform")
{
    xvar = "gmean"
    yvar = "residual_variance"
}

xxm <- melt(xx[, c(xvar,yvar,"var.gene")],
            id.vars=c("var.gene",xvar))

gp <- ggplot(xxm, aes_string(xvar, "value", color="var.gene"))


gp <- gp + scale_color_manual(values=c("black","red"))
gp <- gp + geom_point(alpha = 1, size=0.5)
gp <- gp + facet_wrap(~variable, scales="free")
gp <- gp + theme_classic()

if(opt$normalizationmethod=="sctransform")
{
    gp <- gp + scale_y_continuous(trans="log10", labels = function(x) format(x, digits=3, scientific = TRUE))
    gp <- gp + scale_x_continuous(trans="log10", labels = function(x) format(x, digits=3, scientific = TRUE))
}

gp <- gp + ylab(yvar)

if(exists("xthreshold"))
{
    if(xthreshold > 0)
    {
        gp <- gp + geom_vline(xintercept=xthreshold, linetype="dashed", color="blue")
    }
}

save_ggplots(file.path(opt$outdir, "varGenesPlot"),
             gp=gp,
             width=8,
             height=5)

cat("no. variable genes: ", length(VariableFeatures(object = s)), "\n")
# stats$no_variable_genes <- length(VariableFeatures(object = s))


## ######################################################################### ##
## ###### (ii) Removal of unwanted variation/cell cycle correction ########## ##
## ######################################################################### ##


## initialise the text snippet
tex = ""

## start building figure latex...
subsectionTitle <- getSubsectionTex("Visualisation of cell cycle effects")
tex <- c(tex, subsectionTitle)


## If cell cycle genes are given, make PCAs of the cells based
## on expression of cell cycle genes pre and post normalisation.

cc_genes <- getCellCycleGenes(sgenes_file = opt$sgenes,
                              g2mgenes_file = opt$g2mgenes,
                              species = opt$species)

## score the cell cycle phases
s <- CellCycleScoring(object=s,
                      s.features=cc_genes$s.genes,
                      g2m.features=cc_genes$g2m.genes,
                      set.ident=TRUE)

s <- RunPCA(object = s,
            features = c(cc_genes$s.genes,
                         cc_genes$g2m.genes),
            do.print = FALSE)

## PCA plot on cell cycle genes without regression
gp <- PCAPlot(object = s)

cc_plot_fn <- "cellcycle.without.regression.pca"
cc_plot_path <- file.path(opt$outdir, cc_plot_fn)

save_ggplots(cc_plot_path, gp, width=7, height=4)

pcaCaption <- paste0("PCA analysis of cells based on expression of cell cycle genes ",
                     "(without regression of cell-cyle effects)")

tex <- c(tex, getFigureTex(cc_plot_fn,
                           pcaCaption,
                           plot_dir_var=opt$plotdirvar))

if ( identical(opt$cellcycle, "all") ){

    message("Cell cycle correction for S and G2M scores will be applied\n")
    vars.to.regress <- c(latent.vars, "S.Score", "G2M.Score")

} else if ( identical(opt$cellcycle, "difference") ) {

    message("Cell cycle correction for the difference between G2M and S phase scores will be applied\n")
    s$CC.Difference <- s$S.Score - s$G2M.Score
    vars.to.regress <- c(latent.vars, "CC.Difference")

} else if(opt$cellcycle != "none") {
    stop("Cell cycle regression type not recognised")
}

if(opt$cellcycle != "none")
{
    ## Apply the cell cycle correction.
    ## Always scale the RNA slot

    if(opt$normalizationmethod=="log-normalisation")
    {
        s <- ScaleData(object=s,
                       features = feature.genes,
                       vars.to.regress=vars.to.regress,
                       model.use=opt$modeluse,
                       assay="RNA")

    } else if(opt$normalizationmethod=="sctransform")
    {
        s <- SCTransform(object=s,
                         assay="RNA",
                         new.assay.name="SCT",
                         do.correct.umi=TRUE,
                         variable.features.n=3000,
                         vars.to.regress=vars.to.regress,
                         do.scale=FALSE,
                         do.center=TRUE,
                         return.only.var.genes=FALSE)

        ## Note that the SCT assay will now be the default.
    } else {
        stop("normalisation method not recognised")
    }

    ## visualise the cells by PCA of cell cycle genes after regression
    s <- RunPCA(object = s, features = c(cc_genes$s.genes, cc_genes$g2m.genes), do.print = FALSE)
    gp <- PCAPlot(object = s)

    cc_plot_fn <- paste("cellcycle.regressed", opt$cellcycle, "pca", sep=".")
    cc_plot_path <- file.path(opt$outdir, cc_plot_fn)

    save_ggplots(cc_plot_path, gp, width=7, height=4)

    pcaCaption <- paste0("PCA analysis of cells based on expression of cell cycle genes ",
                         "after regression of cell-cyle effects ",
                         "(regression type: ", opt$cellcycle, ")")

    tex <- c(tex, getFigureTex(cc_plot_fn,
                           pcaCaption,
                           plot_dir_var=opt$plotdirvar))
}

tex_file <- file.path(opt$outdir, "cell.cycle.tex")

writeTex(tex_file, tex)

message("Completed")
