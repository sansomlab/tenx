
# Introduction -----

timestamp()
message("Started")

# Packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
    require(DropletUtils)
    require(rtracklayer)
    require(Matrix)
}))

# Parse options ----

option_list <- list(
    make_option(
        c("--matrixpath", "-m"), action = "store",
        type = "character",
        dest = "matrixpath",
        help="Description of input option"),
    make_option(
        c("--gtf", "-g"), action = "store",
        type = "character",
        dest = "gtf",
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

cat("Loading gene information from GTF file ...")
gtfGeneInfo <- import.gff(opt$gtf, feature.type="gene")
names(gtfGeneInfo) <- gtfGeneInfo$gene_id
cat("Done.\n")

# Warning: the following requires that all "rownames" in the 10x matrix
# are present in the "gene_id" field of the GTF file.
cat("Matching GTF information to 10x matrix ...")
rowRanges(sce10x) <- gtfGeneInfo[rownames(sce10x)]
cat("Done.\n")

cat("Computing total UMIs ...")
sce10x$umis <- colSums(assay(sce10x, "counts"))
cat("Done.\n")

cat("Computing cell ranks ...")
sce10x$rank <- rank(-sce10x$umis, ties.method = "first")
cat("Done.\n")

cat("Computing mitochondrial UMIs ...")
sce10x$umis_mt <- colSums(assay(sce10x[seqnames(sce10x) == "MT", ], "counts"))
cat("Done.\n")

cat("Trimming sample name ...")
sce10x$sample <- gsub("^([[:alnum:]_]+)-count/.*", "\\1", sce10x$Sample)
cat("Done.\n")

cat("Exporting data table ...")
tableExport <- colData(sce10x)[, c("sample", "umis", "rank", "umis_mt")]
write.table(tableExport, opt$outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Done.\n")

# Conclusion ---

message("Completed")
timestamp()
