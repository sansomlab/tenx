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
    make_option(c("--usesigcomponents"), default=FALSE,
                help="use significant principle component"),
    make_option(c("--components"), type="integer", default=10,
                help="if usesigcomponents is FALSE, the number of principle components to use"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use in construction of SNN graph. (e.g. 'pca', 'ica')"),
    make_option(c("--nneighbors"), type="integer", default=20L,
                help="number of neighbors (k.param)"),
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

if(opt$usesigcomponents)
{
    comps <- getSigPC(s)
    message("using the following pcas:")
    print(comps)
} else {
    comps <- 1:as.numeric(opt$components)
}



if(toupper(opt$reductiontype)=="PCA")
{ if(nrow(s@reductions$pca@jackstraw@overall.p.values) > 0)
  {
    ## Make a table of the retained principle components
    x <- as.data.frame(s@reductions$pca@jackstraw@overall.p.values)
    x$p.adj <- p.adjust(x$Score, method="BH")
    x$significant <- "no"
    x$significant[x$p.adj < 0.05] <- "yes"
    x <- x[x$PC %in% comps,]
    x$sdev <- s@reductions$pca@stdev[x$PC]

    print(
        xtable(sprintfResults(x), caption=paste("Table of the selected (n=",
                                                nrow(x),
                                                ") principle components",
                                                sep="")),
        file=file.path(opt$outdir, "selected.principal.components.tex")
    )
  }
}

message("Making the graphs")

#

s <- FindNeighbors(s,
                   reduction = opt$reductiontype,
                   k.param = opt$nneighbors,
                   dims = comps)


message("Finished making the graphs")

message("seurat_begin.R object final default assay: ", DefaultAssay(s))

# Save the graphs
saveRDS(s@graphs, file=file.path(opt$outdir, "graphs.rds"))

message("Completed")
