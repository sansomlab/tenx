## Perform quality assessment of a single-cell RNA-seq dataset 
##
## This script can be compiled from the command line or run interactively in
## Rstudio.


# Libraries --------------------------------------------------------------------
stopifnot(require(yaml))
stopifnot(require(Seurat))
stopifnot(require(sctransform))
stopifnot(require(ggplot2))
stopifnot(require(dplyr))
stopifnot(require(Matrix))
stopifnot(require(xtable))
stopifnot(require(tenxutils))
stopifnot(require(reshape2))
stopifnot(require(future))

# Parameters -------------------------------------------------------------------
#
# The script expects the following paramters:

default_options <- list(
  # Location of the input 10x matrix
  "tenxdir" = "",

  # location for the ouput files
  "outdir" = "",      

  # Path to file containing the cell metadata ("barcode" column must be present) 
  "metadata" = "",       
  
  # The name of a column in the metadata table by which to group samples 
  # (e.g. sample_id)                   # 
  "groupby" = "",       

  "mingenes" = 200,
  
  # Include genes with detected expression in at least this many cells.
  # See Seurat::CreateSeuratObject(min.cells=...)."
  "mincells" = 3,     

  # Whether to randomly downsample cells so that each group has the 
  # same number of cells."
  "downsamplecells" = FALSE, 
  
  # The seed (an integer) to use when down sampling cells
  "seed" = 42, 

  # Minimal count of genes detected to retain a cell. 
  # See Seurat::FilterCells(subset.names='nGene', low.thresholds= ...)
  "qcmingenes" = 500,
 
  # Minimum percentage of UMI assigned to mitochondrial genes. 
  # See Seurat::FilterCells(subset.names='percent.mito', low.thresholds= ...)
  "qcminpercentmito" = -Inf,
  
  # Maximal percentage of UMI assigned to mitochondrial genes.", 
  # See Seurat::FilterCells(subset.names = 'percent.mito', high.thresholds= ...)
  "qcmaxpercentmito" = 0.05, 
  
  # Max nCount_RNA to retain a cell
  "qcmaxcount" = NULL,
  
  # The normalization method to use
  "normalizationmethod" = "log-normalization",

  # Latent variables to regress out.,
  # See Seurat::ScaleData(vars.to.regress=..., model.use=opt$modeluse)
  "latentvars" = "nUMI,percent.mito",

  # Model used to regress out latent variables.,
  # See Seurat::ScaleData(model.use=opt$modeluse)
  "modeluse" = "linear",

  # Method for variable gene selection.
  # Either "top.genes" or "mean.var.plot"         
  "vargenesmethod" = "mean.var.plot",
      
  # Number of highly variable genes to retain      
  "topgenes" = 1000,
  
  # Bottom cutoff on y-axis for identifying variable genes.
  # See Seurat::FindVariableFeatures(y.cutoff=...)
  "sdcutoff" = 0.5,
        help=paste(
          
  # Bottom cutoff on x-axis for identifying variable genes
  # See Seurat::FindVariableFeatures(x.low.cutoff=...)            
  "xlowcutoff" = 0.1,

  # Top cutoff on x-axis for identifying variable genes
  # See Seurat::FindVariableFeatures(x.low.cutoff=...)
  "xhighcutoff" = 8,

  # minimum mean of log counts when using trendvar method
  "minmean" = 0,

  # significance threshold for trendvar method
  "vargenespadjust" = 0.05,
  
  # Whether or not to subset the cells. Either "use.all" or the path of a file 
  # containing the list of barcode ids to retain (no header, 1 per line).
  "subsetcells" = "use.all",

  # A factor specified in metadata.tsv on which to subset
  "subsetfactor" = NULL,

  # The desired level of the sub-setting factor
  "subsetlevel" ="none",

  # A file containing a list of barcode ids to remove (if present)
  # (no header, 1 per line).
  "blacklist" = NULL,

  # type of cell cycle regression to apply (none|all|difference)           
  "cellcycle" = "none",
    
  # A vector of Ensembl gene ids associated with S phases.
  # See Seurat::CellCycleScoring(s.genes=...)
  "sgenes" = "none",
      
  # A vector of Ensembl gene ids associated with G2M phase.
  # See Seurat::CellCycleScoring(g2m.genes=...)
  "g2mgenes" = "none",

  # Number of replicates for the jackstraw analysis     
  "jackstrawnumreplicates" = 100,

  # Number of cores to be used for the Jackstraw analysis
  "numcores" = 12,
        
  # use significant principle component
  "usesigcomponents" = FALSE,

  # if usesigcomponents is FALSE, the number of principle components to use
  "components" = 10
)
 
# If running interactively comment this line out or specify the
# location of a file containing the parameters that you wish to use, e.g.
# options <- read_yaml("quality_assessment.myconfig.yaml")
options <- read_yaml(params@options)

if(! all(names(options) %in% names(default_options))) {
  stop("option not recognised")
}

# Update the default options
if(!is.null(options)) {
  opt <- utils::modifyList(default_options, options)
} else{
  opt <- default_options
}

cat("Running with options:\n")
print(opt)

# Set up for multiprocessing
plan("multiprocess",
     workers = opt$numcores)


# Input data ----

## ######################################################################### ##
## ###################### (i) Read in data ################################# ##
## ######################################################################### ##









## Normalisation specific options.

if(opt$normalizationmethod=="log-normalization")
{

    ## variable gene identification

    message("Finding variable features")
    ## We need to run FindVariableFeatures to set HVFInfo(object = s)
    ## even if "trendvar" method is specified...
    if(opt$vargenesmethod=="trendvar")
    {
        fvg_method="mean.var.plot"
    } else {
        fvg_method=opt$vargenesmethod
    }

    s <- FindVariableFeatures(s,
                              selection.method=fvg_method,
                              nfeatures=opt$topgenes,
                              mean.cutoff = c(opt$xlowcutoff,
                                              opt$xhighcutoff),
                              dispersion.cutoff=c(opt$sdcutoff, Inf))

    xthreshold <- opt$xlowcutoff

    if(opt$vargenesmethod=="trendvar")
    {
        message("setting variable genes using the trendvar method")

        ## get highly variable genes using the getHVG function in tenxutils (Matrix.R)
        ## that wraps the trendVar method from scran.
        hvg.out <- getHVG(s,
                          min_mean=opt$minmean,
                          p_adjust_threshold=opt$vargenespadjust)

        ## overwrite the slot
        VariableFeatures(object = s) <- row.names(hvg.out)

        xthreshold <- opt$minmean
    }


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
    xvar = "mean"
    yvar = "dispersion"

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
stats$no_variable_genes <- length(VariableFeatures(object = s))

print(stats)

# write out some statistics into a latex table.
print(
    xtable(t(data.frame(stats)), caption="Run statistics"),
    file=file.path(opt$outdir, "stats.tex")
    )


## ######################################################################### ##
## ###### (v) Removal of unwanted variation/cell cycle correction ########## ##
## ######################################################################### ##


## initialise the text snippet
tex = ""

## start building figure latex...
subsectionTitle <- getSubsectionTex("Visualisation of cell cycle effects")
tex <- c(tex, subsectionTitle)


## If cell cycle genes are given, make a PCA of the cells based
## on expression of cell cycle genes

cell_cycle_genes <- FALSE

if (!(is.null(opt$sgenes) | opt$sgenes=="none")
    & !(is.null(opt$g2mgenes) | opt$g2mgenes=="none"))
{

  cell_cycle_genes <- TRUE
  cat("There are cell cycle genes")

  # get the genes representing the cell cycle phases
  sgenes_ensembl <- read.table(opt$sgenes, header=F, as.is=T)$V1
  sgenes <- s@misc$seurat_id[s@misc$gene_id %in% sgenes_ensembl]

  g2mgenes_ensembl <- read.table(opt$g2mgenes, header=F, as.is=T)$V1
  g2mgenes <- s@misc$seurat_id[s@misc$gene_id %in% g2mgenes_ensembl]

  # score the cell cycle phases
    s <- CellCycleScoring(object=s,
                          s.features=sgenes,
                          g2m.features=g2mgenes,
                          set.ident=TRUE)

    s <- RunPCA(object = s,
                pc.genes = c(sgenes, g2mgenes),
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

        message("Cell cycle correction for S and G2M scores will be applied\n")
        vars.to.regress <- c(latent.vars, "S.Score", "G2M.Score")

    } else if ( identical(opt$cellcycle, "difference") ) {

      message("Cell cycle correction for the difference between G2M and S phase scores will be applied\n")
        s$CC.Difference <- s$S.Score - s$G2M.Score
        vars.to.regress <- c(latent.vars, "CC.Difference")

    } else {
      stop("Cell cycle regression type not recognised")
    }

    ## Apply the cell cycle correction.

    ## Always scale the RNA slot
    s <- ScaleData(object=s,
                   features = all.genes,
                   vars.to.regress=vars.to.regress,
                   model.use=opt$modeluse,
                   assay="RNA")

    ## Optionally run sctransorm
    if(opt$normalizationmethod=="sctransform")
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
                               plot_dir_var=opt$plotdirvar))

}


tex_file <- file.path(opt$outdir, "cell.cycle.tex")

writeTex(tex_file, tex)


## ######################################################################### ##
## ################### (vi) Dimension reduction (PCA) ###################### ##
## ######################################################################### ##

# perform PCA using the variable genes
s <- RunPCA(s,
            pc.genes=VariableFeatures(object = s),
            npcs=100,
            do.print=FALSE)

n_cells_pca <- min(1000, length(colnames(x = s)))

# Write out heatmaps showing the genes loading the first 12 components
plot_fn <- function() {
    DimHeatmap(s, dims=1:12,
               # cells.use=n_cells_pca,
               reduction="pca",
               balanced=TRUE,
               # label.columns=FALSE,
               # cexRow=0.8,
               # use.full=FALSE
        )
}

save_plots(file.path(opt$outdir, "pcaComponents"), plot_fn=plot_fn,
           width=8, height=12)

print("----")
print(dim(Embeddings(object=s, reduction="pca")))

nPCs <- min(dim(Embeddings(object =s, reduction="pca"))[2],50)

print(nPCs)

# Write out the PCA elbow (scree) plot
png(file.path(opt$outdir, "pcaElbow.png"),
    width=5, height=4, units="in",
    res=300)

ElbowPlot(s,
          ndims=nPCs)

dev.off()

## In Macosko et al, we implemented a resampling test inspired by the jackStraw procedure.
## We randomly permute a subset of the data (1% by default) and rerun PCA,
## constructing a 'null distribution' of gene scores, and repeat this procedure. We identify
## 'significant' PCs as those who have a strong enrichment of low p-value genes.

s <- JackStraw(s,
               reduction="pca",
               num.replicate=opt$jackstrawnumreplicates,
               dims = nPCs)

s <- ScoreJackStraw(s, dims = 1:nPCs)

#               do.par=TRUE,
#               num.cores=opt$numcores)

 gp <- JackStrawPlot(s, dims= 1:nPCs,
                     reduction="pca")

 save_ggplots(paste0(opt$outdir,"/pcaJackStraw"),
            gp,
            width=8,
            height=12)

message("seurat_begin.R object final default assay: ", DefaultAssay(s))

# Save the R object
saveRDS(s, file=file.path(opt$outdir, "begin.rds"))

message("Completed")
