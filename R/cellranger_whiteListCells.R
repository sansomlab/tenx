
# Introduction -----

timestamp()
message("Started")

# Packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
    require(RSQLite)
    require(stringr)
    require(UpSetR)
}))

# Test options ----

opt <- list(
    tablename="cellranger_cellcalling",
    methods="cellranger,emptyDrops",
    combine="union",
    outfile="test.txt.gz"
)

# Parse options ----

option_list <- list(
    make_option(
        c("--tablename", "-t"), action="store",
        type="character",
        dest="tablename",
        help="Name of the table to load from the project database."),
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
        c("--combine", "-c"), action="store",
        type="character",
        dest="combine",
        default="union",
        help=paste(
            "Strategy to combine calls from multiple methods.",
            "Choices are: union, intersect"
        )),
    make_option(
        c("--outfile", "-o"), action="store",
        type="character",
        dest="outfile",
        help="Output PDF file.")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")
print(opt)

# Process data ----

stopifnot(opt$combine %in% c("union", "intersect"))

cat("Identify methods to use ... ")
use_methods <- str_split(string=opt$methods, pattern=",")[[1]]
cat("Done.\n")

cat("Connecting to database ...")
conn <- dbConnect(SQLite(), "csvdb")
cat("Done.\n")

cat("Loading database table ...")
# RSQLite::dbListTables(conn)
tableData <- dbReadTable(conn, opt$tablename)
cat("Done.\n")

dbDisconnect(conn)

cat("Computing white list ...")
timesCalled <- rowSums(tableData[, use_methods])
# table(timesCalled)
combineOperator <- if (opt$combine == "union") {
    " | "
} else if ( opt$combine == "intersect") {
    " & "
} else {
    stop("Invalid `combine` option.")
}
filterExpression <- as.expression(paste(use_methods, collapse = combineOperator))
whiteList <- subset(tableData, eval(parse(text = filterExpression)), c("sample", "barcode"))
cat("Done.\n")

cat("Writing white list ...")
gzfileOut <- gzfile(opt$outfile, open = "wt")
write.table(x=whiteList, file=gzfileOut, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
close(gzfileOut)
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
