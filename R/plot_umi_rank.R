
# Introduction -----

timestamp()
message("Started")

# packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
    require(RSQLite)
    require(ggplot2)
    require(cowplot)
}))

# test options ----

opt <- list(
  tablename = "cellranger_raw_count",
  outfile = "cellranger_raw_count.pdf"
)

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
tableData <- subset(tableData, umis > 0)
cat("Done.\n")

cat("Plotting...")
gg <- ggplot(tableData, aes(rank, umis, color=sample)) +
    geom_line(size=0.2) +
    scale_x_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
    scale_y_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
    ylab("Total UMIs") +
    xlab("Barcode rank") +
    theme(panel.grid.major = element_line(size=0.2, colour="grey"))
ggsave(opt$outfile, gg, width=7, height=4)
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
