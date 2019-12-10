## Title ----
##
## Export the spliced, unspliced and ambiguous matrices ready for
## velocity analysis (e.g. scvelo)
##
## Description ----
##
## The script expects as input a single mm matrix in which genes names
## are prefixed with ex_ (exons), in_ (introns), sp_ (spanning)
##
## The matrices returned have the same dimensions and indices
## (missing values are populated with zeros).
##
## Usage ----
##
## statement = '''Rscript %(tenx_dir)s/R/dropest_export_matrices.R
##                --tenxdir=all-datasets/all.dir
##                --outdir=%(out_dir)s
##                &> %(outfile)s
##             '''

message("dropest_export_matrices.R")
timestamp()

# Libraries ----

stopifnot(
    require(optparse),
    require(R.utils),
    require(Matrix)
)

# Options ----

option_list <- list(
    make_option(
        c("--tenxdir"),
        help="10x matrix directory"
    ),
    make_option(
        c("--outdir"),
        default=".",
        help="location where the outputs will be written"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Functions ----


# Setup output directory ----

if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir)
}


# Sanity checks ----

## Function definitions

## A function to prepare matching matrices ("layers") of spliced, unspliced
## and ambiguous reads for velocity analysis
getLayer <- function(data,
                     feature_names,
                     prefixes,
                     prefix,
                     unique_features)
{
    ## subset layer and features
    loc <- prefixes == prefix
    layer <- data[loc,]
    layer_features = feature_names[loc]

    ## identify features missing from this layer
    to_add = unique_feature_names[!unique_feature_names %in% layer_features]

    missing <- Matrix(0, nrow = length(to_add),
                       ncol = ncol(layer), sparse = TRUE)

    message("adding missing features")
    layer <- rbind(layer, missing)

    rownames(layer) <- c(layer_features,
                         to_add)

    ## reorder the rows to ensure consistency between layers
    layer <- layer[unique_features,]

    layer
}


## Input data ----

## read in the matrix
matrix_path = file.path(opt$tenxdir, "matrix.mtx.gz")
data <- readMM(gzfile(matrix_path))

## read in the features (these are expected to be prefixed with ex_, in_, sp_)
features_path <- file.path(opt$tenxdir, "features.tsv.gz")
features <- read.table(gzfile(features_path), header=F, sep="\t", as.is=T)$V1
prefixes <- substring(features,1,3)
feature_names <- substring(features,4)
unique_feature_names <- unique(feature_names)

## read in the metadata
metadata_path <- file.path(opt$tenxdir, "metadata.tsv.gz")
stopifnot(file.exists(metadata_path))

metadata <- read.table(gzfile(metadata_path), header=TRUE)

## prepare the separate matrices ("layers") for exons, introns and
## spanning reads
layers <- list(exons="ex_",
               introns="in_",
               spanning="sp_")

for(layer in names(layers))
{
    message("exporting layer:", layer)

    xx <- getLayer(data, feature_names,
                    prefixes, layers[layer],
                    unique_feature_names)

    message("layer dimensions:")
    print(dim(xx))

    message("saving layer: ", layer)
    ## write out the layer
    matrix_path = file.path(opt$outdir, paste(layer,"mtx", sep="."))
    writeMM(xx, matrix_path)
    gzip(matrix_path, overwrite=TRUE)
}

## write out the feature names
feature_path = file.path(opt$outdir, "features.tsv")

features = data.frame(feature=unique_feature_names)
print(head(features))
write.table(features, feature_path,
            row.names=F, col.names=F, quote=F, sep="\t")
gzip(feature_path, overwrite=TRUE)

## copy over the cell barcodes
barcodes_path = file.path(opt$outdir, "barcodes.tsv.gz")
file.copy(file.path(opt$tenxdir, "barcodes.tsv.gz"),
          barcodes_path)

## copy over the metadata
metadata_out_path = file.path(opt$outdir, basename(metadata_path))
file.copy(metadata_path, metadata_out_path)

timestamp()
message("Completed")
