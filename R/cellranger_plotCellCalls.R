
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
    outfile="test.pdf"
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
        c("--outfile", "-o"), action="store",
        type="character",
        dest="outfile",
        help="Output PDF file.")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")
print(opt)

# Process data ----

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

cat("Computing sets ...")
whichInSet <- function(i, x){
    which(x[[i]] == 1)
}
upsetDataList <- lapply(use_methods, whichInSet, x=tableData)
names(upsetDataList) <- use_methods
cat("Done.\n")

cat("Plotting ...")
pdf(opt$outfile, width = 7, height = 5)
upset(data = fromList(upsetDataList), order.by = "freq")
dev.off()
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
