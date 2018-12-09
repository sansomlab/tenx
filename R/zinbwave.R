## Run zinbwave on a seurat object.

# Libraries ----

stopifnot(
  require(Seurat),
  require(tenxutils),
  require(optparse),
  require(zinbwave),
  require(BiocParallel)
)

# Options ----

## deal with the options
option_list <- list(
    make_option(c("--seuratobject"), default="must_specify",
                help="location of the seurat object"),
    make_option(c("--K"), default=10,
                help="zinbwave K parameter"),
    make_option(c("--X"), default="none",
                help="zinbwave X parameter, must match column in metadata"),
    make_option(c("--ncpu"), default=8,
                help="ncpu to use"),
    make_option(c("--ngenes"), type="integer", default=-1,
                help="use a random sample of n genes"),
    make_option(c("--ncells"), type="integer", default=-1,
                help="use a random sample of n cells to use"),
    make_option(c("--backend"),default="snow",
                help="multi,dopar,snow,serial"),
    make_option(c("--outfile"),default="none",
                help="The name of a file to save the modified seurat object to")
    )

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with options:")
print(opt)

## need to module unload apps/gsl!
s <- readRDS(opt$seuratobject)
data <- s@raw.data[s@var.genes,s@cell.names]

message("size of the data matrix:")
print(dim(data))

## random downsampling of the matrix is supported for testing
## purposes..
if(opt$ngenes != -1)
{
    message("subsetting genes")
    genes.take <- sample(rownames(data), opt$ngenes)
    data <- data[genes.take, ]
    print(dim(data))
}
if(opt$ncells != -1)
{
    message("subsetting cells")
    ## draw random sample to try and include all X levels..
    cells.take <- sample(colnames(data), opt$ncells)
    data <- data[, cells.take]
    print(dim(data))
}

# see: https://github.com/drisso/zinbwave/issues/17
se <- SummarizedExperiment(as.matrix(data),colData=s@meta.data[colnames(data),])

print(colnames(colData(se)))

if(opt$backend == "serial")
{
    engine <- SerialParam()
} else if(opt$backend == "multi") {
    engine <- MulticoreParam(workers = opt$ncpu)
} else if(opt$backend == "dopar") {
    engine <- DoparParam(workers = opt$ncpu)
} else if(opt$backend == "snow") {
    engine <- SnowParam(workers = opt$ncpu, type = "SOCK")
} else {
    stop("backend not recognised")
}

message("running zinbwave")
if(opt$X!="none")
{
    X=paste0("~",opt$X)
    print(paste("X specified as:",X))
    z <- zinbwave(se, X=as.formula(X),K=opt$K, BPPARAM=engine)
} else {
    z <- zinbwave(se, K=opt$K, BPPARAM=engine)
}

message("adding zinbwave information to seurat object")
W <- reducedDim(z)

print(head(W))

if(opt$outfile != "none")
{

    if(opt$ncells != -1)
    {
        s <- SubsetData(s,
                        cells.use=cells.take)
    }

    ## set zinbwave as a dim reduction in Seurat
    s <- SetDimReduction(object = s,
                         reduction.type = "zinbwave",
                         slot = "cell.embeddings",
                         new.data = W)

    s <- SetDimReduction(object = s,
                         reduction.type = "zinbwave",
                         slot = "key",
                         new.data = "zinbwave")


    message("saving the R object")

    ## Save the R object
    saveRDS(s,file=opt$outfile)
}

message("all done")
