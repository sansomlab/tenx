
# Introduction -----

timestamp()
message("Started")

# packages ----

stopifnot(suppressPackageStartupMessages({
  require(optparse)
  require(DropletUtils)
  require(Matrix)
}))

# test options ----

opt <- list(
  matrixpath = "donor1_butyrate-count/outs/raw_gene_bc_matrices/GRCh38",
  outfile = "donor1_butyrate-count/cellranger.raw.counts.txt"
)

# Parse options ----

option_list <- list(
  make_option(
    c("--matrixpath", "-m"), action = "store",
    type = "character",
    dest = "matrixpath",
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

cat("Loading matrix data ...")
sce10x <- read10xCounts(opt$matrixpath)
cat("Done.\n")

cat("Computing total UMIs ...")
sce10x$umis <- colSums(assay(sce10x, "counts"))
cat("Done.\n")

cat("Computing cell ranks ...")
sce10x$rank <- rank(-sce10x$umis, ties.method = "first")
cat("Done.\n")

cat("Trimming sample name ...")
sce10x$sample <- gsub("^([[:alnum:]_]+)-count/.*", "\\1", sce10x$Sample)
cat("Done.\n")

cat("Exporting data table ...")
tableExport <- colData(sce10x)[, c("sample", "umis", "rank")]
write.table(tableExport, opt$outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Done.\n")

# Test code for plotting (move to other task)

if (interactive()) {
    ggData <- colData(sce10x)[, c("umis", "rank")]
    ggData <- as.data.frame(ggData)
    ggData <- subset(ggData, umis > 0)
    require(ggplot2)
    require(cowplot)
    ggplot(ggData) +
      geom_line(aes(rank, umis), size=0.2) +
      scale_x_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6)) +
      scale_y_log10(breaks=10^seq(0, 6), labels=10^seq(0, 6))
}

# Conclusion ---

message("Completed")
timestamp()
