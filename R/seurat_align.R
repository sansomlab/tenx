## Cluster the single cells in a given Seurat object

# Libraries ----

stopifnot(
  require(Seurat),
  require(tenxutils),
  require(scater),
  require(scran),
  require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.rds",
                help="A single seurat object containing all the samples to be aligned"),
    make_option(c("--metavar"), default=NULL,
                help="The meta data variable that stratifies the samples to be aligned"),
    make_option(c("--strata"), default="all",
                help="the levels of the meta var that should be included in the alignment"),
    make_option(c("--niter"), type="integer", default=50,
                help="number of iternations"),
    make_option(c("--numccs"), type="integer", default=1,
                help="number of ccs to return"),
    make_option(c("--outfile"), default="aligned.rds",
                help="the name for the output file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Functions ----
getVarGenes <- function(s)
{
    ss <- as.SingleCellExperiment(s)
    var.fit.nospike <- trendVar(ss, parametric=TRUE,
                                use.spikes=FALSE, loess.args=list(span=0.2))
    var.out.nospike <- decomposeVar(ss, var.fit.nospike,
                                    subset.row=rowMeans(as.matrix(logcounts(ss))) > 0.01)
    g.out <- var.out.nospike[which(var.out.nospike$FDR <= 0.05),]
    rm(ss)
    row.names(g.out)
}


runMultiCCA <- function(s,
                        metavar="condition",
                        strata=NULL,
                        niter=50,
                        num.ccs=1)
{
    objects <- list()
    var.genes <- c()
    for(stratum in strata)
    {
        message("subsetting to ", stratum)
        objects[[stratum]] <- SubsetData(s,
                                         cells.use=rownames(s@meta.data)[s@meta.data[[metavar]]==stratum],
                                         subset.raw=T)
        print(dim(objects[[stratum]]@data))

        message("getting variable genes")
        vg <- getVarGenes(objects[[stratum]])
        print(length(vg))
        var.genes <- unique(c(var.genes, vg))
    }

    message("number of variable genes (union):")
    print(length(var.genes))

    message("running alignment")
    RunMultiCCA(objects,
                niter=niter,
                genes.use=var.genes, num.ccs=num.ccs)
}

s <- readRDS(opt$seuratobject)

if(opt$strata=="all")
{
    strata <- unique(s@meta.data[[opt$metavar]])
} else {
    strata <- as.vector(unlist(strsplit(opt$strata, ",")))
}

if(!length(strata)>1)
{
    stop("More than one stratum must be specified")
}

message("strata defined as:")
print(strata)

alignment <- runMultiCCA(s, metavar=opt$metavar,
                   strata=strata,
                   niter=opt$niter,
                   num.ccs=opt$numccs)

message("saving result")

saveRDS(alignment,
        file=opt$outfile)

message("Completed")
