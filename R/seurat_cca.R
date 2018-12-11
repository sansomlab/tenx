## Perform the canonical clustering analysis

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
        hvg.out <- getHVG(objects[[stratum]], min_mean=0.01)
        vg <- rownames(hvg.out)

        print(length(vg))
        var.genes <- unique(c(var.genes, vg))
    }

    message("number of variable genes (union):")
    print(length(var.genes))

    if(length(strata)>2)
    {
        message("running multi-cca alignment")
        RunMultiCCA(objects,
                    niter=niter,
                    genes.use=var.genes, num.ccs=num.ccs)
    } else {
        message("running single cca alignment")
        RunCCA(objects[[strata[1]]],
               objects[[strata[2]]],
               genes.use = var.genes,
               num.cc = num.ccs)
    }
}

s <- readRDS(opt$seuratobject)

all_strata <- unique(s@meta.data[[opt$metavar]])

if(opt$strata=="all")
{
    message("running with all strata:")
    strata <- all_strata
    cells.take <- s@cell.names
} else {
    message("only the given strata will be used:")
    strata <- as.vector(unlist(strsplit(opt$strata, ",")))
    cells.take <- s@meta.data$barcode[s@meta.data[[opt$metavar]] %in% strata]
    if(length(strata[strata %in% all_strata]) != length(strata))
    {
        stop("Not all of the given strata were found in the object")
    }
}

print(strata)

if(!length(strata)>1)
{
    stop("More than one stratum must be specified")
}

cca <- runMultiCCA(s, metavar=opt$metavar,
                   strata=strata,
                   niter=opt$niter,
                   num.ccs=opt$numccs)


message("saving result")

saveRDS(cca,
        file=opt$outfile)

message("Completed")
