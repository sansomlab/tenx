
# Introduction -----

timestamp()
message("Started")

# Packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
    require(stringr)
    require(DropletUtils)
}))

# Test options ----

opt <- list(
    matrixpath="donor1_butyrate-count/outs/raw_gene_bc_matrices/GRCh38",
    methods="cellranger,emptyDrops",
    outfile="donor1_butyrate-count/callCells.txt",
    emptydropslower=100,
    emptydropsniters=10000,
    emptydropsambient=FALSE,
    emptydropsignore=NULL,
    emptydropsretain=NULL,
    emptydropsfdr=0.01,
    threads=4
)

# Parse options ----

option_list <- list(
    make_option(
        c("--matrixpath", "-t"), action="store",
        type="character",
        dest="matrixpath",
        help="Path to the input 10x matrix directory."),
    make_option(
        c("--methods", "-m"), action="store",
        type="character",
        dest="methods",
        default="cellranger",
        help=paste(
            "Comma-separated list of methods required.",
            "Choices are: cellranger, emptyDrops."
        )),
    make_option(
        c("--outfile", "-o"), action="store",
        type="character",
        dest="outfile",
        help="Output TAB-separated file."),
    make_option(
        c("--emptydropslower", "-l"), action="store",
        type="integer",
        dest="emptydropslower",
        default=100,
        help="A numeric scalar specifying the lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets."),
    make_option(
        c("--emptydropsniters", "-n"), action="store",
        type="integer",
        dest="emptydropsniters",
        default=10000,
        help="An integer scalar specifying the number of iterations to use for the Monte Carlo p-value calculations."),
    make_option(
        c("--emptydropsambient", "-a"), action="store_true",
        type="logical",
        dest="emptydropsambient",
        default = FALSE,
        help="A logical scalar indicating whether results should be returned for barcodes with totals less than or equal to lower."),
    make_option(
        c("--emptydropsignore", "-i"), action="store",
        type="integer",
        dest="emptydropsignore",
        default=NULL,
        help="A numeric scalar specifying the lower bound on the total UMI count, at or below which barcodes will be ignored (see Details for how this differs from lower)."),
    make_option(
        c("--emptydropsretain", "-r"), action="store",
        type="integer",
        dest="emptydropsretain",
        default=NULL,
        help="A numeric scalar specifying the threshold for the total UMI count above which all barcodes are assumed to contain cells."),
    make_option(
        c("--emptydropsfdr", "-f"), action="store",
        type="double",
        dest="emptydropsfdr",
        default=0.01,
        help="FDR to call cells that differ from ambient profile."),
    make_option(
        c("--threads", "-p"), action="store",
        type="integer",
        dest="threads",
        default=1,
        help="Number of parallel threads to compute p-values.")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")
print(opt)

# Process data ----

cat("Identify methods to use ... ")
use_methods <- str_split(string=opt$methods, pattern=",")[[1]]
cat("Done.\n")

cat("Import 10x count data ... ")
sce <- read10xCounts(opt$matrixpath, col.names=TRUE)
cat("Done.\n")

cat("Initialize output table ... ")
cellCallTables <- data.frame(
    sample=gsub("(.+)-count\\/.*", "\\1", opt$matrixpath),
    barcode=colnames(sce),
    stringsAsFactors=FALSE
)
cat("Done.\n")

currentMethod <- "cellranger"
if (currentMethod %in% use_methods) {
    cat("Import barcodes from CellRanger filtered matrix ... ")
    filteredMatrixFolder <- gsub("\\/raw_gene_bc_matrices\\/", "/filtered_gene_bc_matrices/", opt$matrixpath)
    stopifnot(dir.exists(filteredMatrixFolder))
    # list.files(filteredMatrixFolder)
    infile <- file.path(filteredMatrixFolder, "barcodes.tsv")
    cellRangerCellBarcodes <- scan(file=infile, what="character")
    cellCallTables[, currentMethod] <- cellCallTables$barcode %in% cellRangerCellBarcodes
    cat("Done.\n")
}

currentMethod <- "emptyDrops"
if (currentMethod %in% use_methods) {
    cat(sprintf("Call cells using `%s` ... ", currentMethod))
    bp <- SerialParam()
    if (opt$threads > 1) {
        bp <- MulticoreParam(workers=opt$threads)
    }
    out <- emptyDrops(
        m=assay(sce, "counts"), lower=opt$emptydropslower, retain=opt$emptydropsretain,
        niters=opt$emptydropsniters, test.ambient=opt$emptydropsambient, ignore=opt$emptydropsignore,
        BPPARAM=bp)
    # table(out$FDR <= opt$emptydropsfdr)
    cellCallTables[, currentMethod] <- (out$FDR <= opt$emptydropsfdr)
    cat("Done.\n")
}
# with(cellCallTables, table(cellranger, emptyDrops, useNA="ifany"))

cat(sprintf("Write table of cell calls ... ", currentMethod))
gzfileOut <- gzfile(opt$outfile, open = "wt")
write.table(x=cellCallTables, file=gzfileOut, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
close(gzfileOut)
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
