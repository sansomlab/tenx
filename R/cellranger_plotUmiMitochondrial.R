
# Introduction -----

timestamp()
message("Started")

# Packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
    require(RSQLite)
    require(ggplot2)
    require(cowplot)
    require(scales)
}))

# Parse options ----

option_list <- list(
    make_option(
        c("--tablename", "-m"), action = "store",
        type = "character",
        dest = "tablename",
        help="Description of input option"),
    make_option(
        c("--outfile", "-o"), action = "store",
        type = "character",
        dest = "outfile",
        help="Description of input option")
)

opt <- parse_args(OptionParser(option_list = option_list))

message("Running with options:")
print(opt)

# Process data ----

cat("Connecting to database ...")
conn <- RSQLite::dbConnect(SQLite(), "csvdb")
cat("Done.\n")

cat("Loading database table ...")
# RSQLite::dbListTables(conn)
tableData <- dbReadTable(conn, opt$tablename)
cat("Done.\n")

dbDisconnect(conn)

cat("Preprocessing data for ggplot ...")
tableData <- subset(tableData, umis > 10)
tableData$umis_mt_percent <- tableData$umis_mt / tableData$umis
cat("Done.\n")

cat("Plotting...")
numberSamples <- length(unique(tableData$sample))
gg <- ggplot(tableData, aes(umis_mt_percent, umis)) +
    facet_wrap(~sample, ncol=2) +
    stat_binhex(aes(fill=log10(..count..))) +
    scale_x_continuous(labels = percent_format()) +
    scale_y_log10(breaks=10^seq(1, 6), labels=10^seq(1, 6)) +
    scale_fill_viridis_c(labels=function(x){ 10^x }) +
    ylab("Total UMIs") +
    xlab("Fraction of mitochondrial UMIs") +
    labs(fill=expression("Barcodes ("*log[10]*")")) +
    theme(panel.grid.major = element_line(size=0.2, colour="grey"))
ggsave(opt$outfile, gg, width=7, height=2*ceiling(numberSamples/2))
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
