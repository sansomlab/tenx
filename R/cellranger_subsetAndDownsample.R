## Title ----
##
## Subset the aggregated and cleaned matrix to subsets of samples
##
## Description ----
##
## A script to downsample the clean aggregated matrix generated with the
## cellranger_subset script.
##
## The script can downsample UMI counts across a (sub)set of samples
## to match their median (for instance).
##
## Usage ----
##
## statement = '''Rscript %(tenx_dir)s/R/downsampleCleanMatrix.R
##                --sampleids=%(sampleIds)s
##                --tenxdir=all-clean/agg.clean.dir
##                --downsample=median
##                --apply=after_subsetting
##                --outdir=%(out_dir)s
##                &> %(outfile)s
##             '''

message("cellranger_subsetAndDownsample.R")
timestamp()

# Libraries ----

stopifnot(
  require(optparse),
  require(methods), # https://github.com/tudo-r/BatchJobs/issues/27
  require(Matrix),
  require(S4Vectors),
  require(ggplot2),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
        c("--tenxdir"),
        help="10x matrix directory"
    ),
    make_option(
        c("--sampleids"),
        help="Comma separated list of sample identifiers in metadata.tsv"
    ),
    make_option(
        c("--downsample"),
        default="median",
        help="Strategy to normalise UMI between multiple 10x samples"
    ),
    make_option(
        c("--apply"),
        default="after_subsetting",
        help=paste('When to apply downsampling. Must be set to either: ',
                   '"after_subsetting" or "before_subsetting"'
        )
    ),
    make_option(
        c("--samplenamefields"),
        default="sample",
        help=paste(
            "The sample name fields supplied as a comma separated list.",
            "Sample name fields must be separated by underscores.",
            "See also script cellranger_cleanAggMatrix.R."
        )
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

downsample <- function(matrixUMI)
{

  message("making pre-downsampling diagnostic plot")
  plotDownsampling(matrixUMI, metadata, "UMI_input")

  agg_ids <- barcode2table(colnames(matrixUMI))$agg_id

  message("downsampling the matrix")
  matrixUMI <- downsampleMatrix(matrixUMI,
                                downsample_method=opt$downsample,
                                library_ids=agg_ids)

  message("making post-downsampling diagnostic plot")
  plotDownsampling(matrixUMI, metadata, "UMI_downsampled")
  matrixUMI
}


# Setup output directory ----

if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir)
}


# Sanity checks ----

if(!opt$apply %in% c("after_subsetting", "before_subsetting"))
{
    stop(paste('--apply must be set to either "after_subsetting"',
               'or "before_subsetting"'))
}


# Input data ----

## Matrix

matrixFile <- file.path(opt$tenxdir, "matrix.mtx.gz")
stopifnot(file.exists(matrixFile))
cat("Importing matrix from:", matrixFile, " ... ")
matrixUMI <- readMM(gzfile(matrixFile))

cat("Done.\n")
cat(
    "Input matrix size:",
    sprintf("%i rows/genes, %i columns/cells\n", nrow(matrixUMI), ncol(matrixUMI))
)

## Barcodes
barcodeFile <- file.path(opt$tenxdir, "barcodes.tsv.gz")
stopifnot(file.exists(barcodeFile))
cat("Importing cell barcodes from:", barcodeFile, " ... ")
barcodes <- scan(gzfile(barcodeFile), "character")

## Metadata

metadataFile <- file.path(opt$tenxdir, "metadata.tsv.gz")
stopifnot(file.exists(metadataFile))
cat("Importing metadata from:", metadataFile, " ... ")
metadata <- read.table(gzfile(metadataFile), header=TRUE)
cat("Done.\n")
cat(
    "Input metadata size:",
    sprintf("%i rows/cells, %i columns\n", nrow(metadata), ncol(metadata)),
    "\n"
)

# Preprocess ----

## Metadata
metadata$code <- as.character(metadata$code)
metadata$agg_id <- as.factor(metadata$agg_id)
metadata$barcode <- as.character(metadata$barcode)
metadata$seq_id <- as.factor(metadata$seq_id)

## UMI matrix colnames / cell barcodes
colnames(matrixUMI) <- barcodes

# Apply downsampling before subsetting, if required ----
if (!identical(opt$downsample, "no") && opt$apply == "before_subsetting") {

  cat("Applying downsampling before subsetting.\n")
  matrixUMI <- downsample(matrixUMI)
}

# Subset ----

sampleIds <- strsplit(opt$sampleids, ",")[[1]]

if(!all(sampleIds %in% metadata$sample_id))
{
    print(paste("given sample ids:", paste(sampleIds, collapse=",")))
    print(paste("metadata sample ids:",
          paste(unique(metadata$sample_id), collapse=",")))
    stop("Not all sampleIds present in the medatadata")
}

if (length(unique(sampleIds)) < length(sampleIds)) {
    print(sampleIds)
    stop("A sample id is present multiple times")
}

if (length(sampleIds) == nlevels(metadata$sample_id)) {
    if (!identical(basename(opt$outdir), "all")) {
        warning(paste(
            "Complete set of samples required.",
            "Basename of --outdir should be 'all'"
        ))
    }
}

# Identify cells to keep
keepCells <- metadata$sample_id %in% sampleIds

# Print information
cat("Count of cells kept (TRUE):\n")
print(table(keepCells))
keepProportion <- sum(keepCells) /length(keepCells)
cat(sprintf("Proportion of cells kept: %.1f%%\n", 100* keepProportion))

# Subset
matrixUMI <- matrixUMI[, keepCells]
barcodes <- barcodes[keepCells]
metadata <- droplevels(metadata[keepCells,])

# Print information
cat(
    "New matrix size (cells subset):",
    sprintf("%i rows/genes, %i columns/cells\n", nrow(matrixUMI), ncol(matrixUMI))
)
cat(sprintf("New count of barcodes: %i", length(barcodes)))
cat(
    "New metadata size:",
    sprintf("%i rows/cells, %i columns\n", nrow(metadata), ncol(metadata))
)

# Apply downsampling after subsetting, if required ----
if (!identical(opt$downsample, "no") && opt$apply == "after_subsetting") {

  cat("Applying downsampling after subsetting.\n")
  matrixUMI <- downsample(matrixUMI)
}

# Write out matrix ----

genesFile <- file.path(opt$tenxdir,"features.tsv.gz")
writeMatrix(opt$outdir, matrixUMI, barcodes, genesFile, metadata)

timestamp()
message("Completed")
