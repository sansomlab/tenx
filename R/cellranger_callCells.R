
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
    tenxdir="donor1_butyrate-count/outs/raw_gene_bc_matrices/GRCh38",
    methods="cellranger",
    outfile="donor1_butyrate-count/callCells.txt"
)

# Parse options ----

option_list <- list(
    make_option(
        c("--tenxdir", "-t"), action="store",
        type="character",
        dest="tenxdir",
        help="Path to the input 10x matrix directory."),
    make_option(
        c("--methods", "-m"), action="store",
        type="character",
        dest="methods",
        default="cellranger",
        help=paste(
            "Comma-separated list of methods required.",
            "Choices are: cellranger."
        )),
    make_option(
        c("--outfile", "-o"), action="store",
        type="character",
        dest="outfile",
        help="Output TAB-separated file."),
    make_option(
        c("--longflag", "-s"), action="store",
        type=c("character", "integer", "logical", "double", "complex"),
        dest="long_flag",
        default="if_applicable",
        help="Description of input option")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")
print(opt)

# Process data ----

cat("Identify methods to use ... ")
use_methods <- str_split(string=opt$methods, pattern=",")[[1]]
cat("Done.\n")

cat("Import 10x count data ... ")
sce <- read10xCounts(opt$tenxdir, col.names=TRUE)
cat("Done.\n")

cat("Initialize output table ... ")
cellCallTables <- data.frame(
    sample=gsub("(.+)-count\\/.*", "\\1", opt$tenxdir),
    barcode=colnames(sce),
    stringsAsFactors=FALSE
)
cat("Done.\n")

currentMethod <- "cellranger"
if (currentMethod %in% use_methods) {
    cat("Import barcodes from CellRanger filtered matrix ... ")
    filteredMatrixFolder <- gsub("\\/raw_gene_bc_matrices\\/", "/filtered_gene_bc_matrices/", opt$tenxdir)
    stopifnot(dir.exists(filteredMatrixFolder))
    # list.files(filteredMatrixFolder)
    infile <- file.path(filteredMatrixFolder, "barcodes.tsv")
    cellRangerCellBarcodes <- scan(file=infile, what="character")
    cellCallTables[, currentMethod] <- cellCallTables$barcode %in% cellRangerCellBarcodes
    cat("Done.\n")
}

write.table(x=cellCallTables, file=opt$outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Conclusion ---

message("Completed")
timestamp()
