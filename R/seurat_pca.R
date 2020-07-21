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
        default="none",
        help=paste("rds file containing the seurat object")
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
        c("--vargenesmethod"),
        default="mean.var.plot",
        help=paste(
            "Method for variable gene selection.",
            "Either top.genes or mean.var.plot"
            )
        ),
    make_option(
        c("--jackstraw"),
        action = "store_true",
        default=FALSE,
        help="should the jackstraw analysis be run"
    ),

    make_option(
        c("--jackstrawnumreplicates"),
        type="integer",
        default=100,
        help="Number of replicates for the jackstraw analysis"
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

s <- loadSeurat(path=opt$seuratobject)

## ######################################################################### ##
## ################### (vi) Dimension reduction (PCA) ###################### ##
## ######################################################################### ##

message("Running the PCA")

# perform PCA using the variable genes
s <- RunPCA(s,
            features=VariableFeatures(object = s),
            npcs=100,
            verbose = FALSE,
            ndims.print = 1:5,
            nfeatures.print = 5)

message("PCA complete")

gc()

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

if(opt$normalizationmethod!="sctransform" & opt$jackstraw)
{

    message("Running JackStraw analysis")

    s <- JackStraw(s,
                   reduction="pca",
                   num.replicate=opt$jackstrawnumreplicates,
                   dims = nPCs)

    s <- ScoreJackStraw(s, dims = 1:nPCs)

    ##               do.par=TRUE,
    ##               num.cores=opt$numcores)

    gp <- JackStrawPlot(s, dims= 1:nPCs,
                     reduction="pca")

    save_ggplots(paste0(opt$outdir,"/pcaJackStraw"),
                 gp,
                 width=8,
                 height=12)
}
message("seurat_begin.R object final default assay: ", DefaultAssay(s))

# Save the R object
saveSeurat(path=opt$seuratobject)

message("Completed")
