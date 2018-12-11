## Align the cca subspace

# Libraries ----

stopifnot(
  require(Seurat),
  require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.rds",
                help="A single seurat object containing all the samples to be aligned"),
    make_option(c("--metavar"), default=NULL,
                help="The meta data variable that stratifies the samples to be aligned"),
    make_option(c("--dimsalign"), type="integer", default=1,
                help="number of ccs to return"),
    make_option(c("--outfile"), default="aligned.rds",
                help="the name for the output file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


cca <- readRDS(opt$seuratobject)

message("Running alignment")

cca.aligned <- AlignSubspace(cca,
                             reduction.type = "cca",
                             grouping.var = opt$metavar,
                             dims.align = 1:opt$dimsalign)


message("saving result")

saveRDS(cca.aligned,
        file=opt$outfile)

message("Completed")
