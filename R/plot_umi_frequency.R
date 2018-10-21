
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
conn <- dbConnect(SQLite(), "csvdb")
cat("Done.\n")

cat("Loading database table ...")
# RSQLite::dbListTables(conn)
tableData <- dbReadTable(conn, opt$tablename)
cat("Done.\n")

dbDisconnect(conn)

cat("Preprocessing data for ggplot ...")
tableData <- subset(tableData, umis > 100)
cat("Done.\n")

cat("Plotting...")
gg <- ggplot(tableData, aes(umis, color=sample)) +
  # geom_density(size=0.2) +
  geom_histogram(
    position = position_identity(),
    alpha=0) +
  facet_wrap(~sample, ncol=1) +
  scale_x_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
  xlab("Total UMIs") +
  ylab("Barcodes") +
  theme(panel.grid.major = element_line(size=0.2, colour="grey"))
ggsave(opt$outfile, gg, width=7, height=2*length(unique(tableData$sample)))
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
