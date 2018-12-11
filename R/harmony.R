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
    make_option(c("--metavar"), default=NULL,
                help="A meta data variable that stratifies the samples"),
    make_option(c("--strata"), default="all",
                help="the levels of the meta var that should be included in the analysis"),
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

message("dimensions of loaded seurat raw data")
print(dim(s@raw.data))


## Step 1: Subset to cells of interest (if necessary)
## --------------------------------------------------

# 1.1 determine the cells that belong to the strata of interest

all_strata <- unique(s@meta.data[[opt$metavar]])

if(opt$strata=="all")
{
    message("Keeping all strata:")
    strata <- all_strata
    cells.take <- s@cell.names
} else {
    message("Using only the given strata:")
    strata <- as.vector(unlist(strsplit(opt$strata, ",")))
    cells.take <- rownames(s@meta.data)[s@meta.data[[opt$metavar]] %in% strata]
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

message("number of cells in the selected strata:")
print(length(cells.take))


## 1.2 if ncells is specified, randomly downsample to the given numbers
## ! note that this does not guarentee that all strata will be represented !

if(opt$ncells != -1 & opt$ncells < length(cells.take))
{
    message("downsampling cell number")
    ## draw random sample to try and include all X levels..
    cells.take <- sample(cells.take, opt$ncells)

    message("number of cells after downsampling:")
    print(length(cells.take))

}


## ensure that we have the right data
if(cells.take != s@cell.names)
{
    s <- SubsetData(s,
                    cells.use=cells.take,
                    subset.raw=TRUE)

}


message("dimensions of seurat raw data after strata selection")
print(dim(s@raw.data))


## 2: Select the genes that will go into the analysis
## --------------------------------------------------


message("Identifying variable genes (trend var) method")
## ensure that variable genes are identified from the
## the data that will be used.
hvg.out <- getHVG(s,
                  min_mean=0,
                  p_adjust_threshold=0.05)

message("number of variable genes:")
print(nrow(hvg.out))
print(head(hvg.out))

genes.take <- rownames(hvg.out)
print(head(genes.take))

## purposes..
if(opt$ngenes != -1 & opt$ngenes < length(genes.take))
{
    ## Here genes are sorted from most to least variable
    ## so we keep the most variable
    message("subsetting genes")
    genes.take <- genes.take[1:opt$ngenes]
}


## 3: Get the raw data for zimbwave
## --------------------------------

## Extract the raw data from the Seurat object
data <- s@raw.data[genes.take,
                   cells.take]

message("size of the data matrix:")
print(dim(data))

## Make the summarised experiement object

# see: https://github.com/drisso/zinbwave/issues/17
se <- SummarizedExperiment(as.matrix(data),
                           colData=s@meta.data[colnames(data),])

print(colnames(colData(se)))

## 4: Run Zinbwave
## ----------------

## Set the parallelisaiton backend
##
## snow has issues running on the cluster so use of
## multicoreparam is now enforced
##
if(opt$backend == "serial")
{
    engine <- SerialParam()
} else if(opt$backend == "multi") {
    engine <- MulticoreParam(workers = opt$ncpu)
} else if(opt$backend == "dopar") {
    engine <- DoparParam(workers = opt$ncpu)
} else {
    stop('only the "serial", "multi" and "dopar" backends are supported')
}

## Run zinbwave
message("running zinbwave")
if(opt$X!="none")
{
    # build the formula for X
    X=paste0("~",opt$X)
    print(paste("X specified as:",X))

    z <- zinbwave(se,
                  X=as.formula(X),
                  K=opt$K,
                  verbose=TRUE,
                  BPPARAM=engine)

} else {
    # run without the X parameter
    z <- zinbwave(se, K=opt$K, BPPARAM=engine)
}

message("adding zinbwave information to seurat object")
W <- reducedDim(z)

print(head(W))

if(opt$outfile != "none")
{

    ## set zinbwave as a dim reduction in Seurat
    s <- SetDimReduction(object = s,
                         reduction.type = "zinbwave",
                         slot = "cell.embeddings",
                         new.data = W)

    s <- SetDimReduction(object = s,
                         reduction.type = "zinbwave",
                         slot = "key",
                         new.data = "zinbwave")

    s@var.genes <- genes.take

    message("saving the R object")

    ## Save the R object
    saveRDS(s,
            file=opt$outfile)
}

message("all done")
