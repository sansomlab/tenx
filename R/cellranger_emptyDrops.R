
# Introduction -----

timestamp()
message("Started")

# Packages ----

stopifnot(suppressPackageStartupMessages({
    require(optparse)
    require(DropletUtils)
}))

# Test options ----

opt <- list(
    tenxdir="donor1_butyrate-count/outs/raw_gene_bc_matrices/GRCh38",
    fdr=0.01,
    outfile="donor1_butyrate-count/DropletUtils_emptyDrops.txt"
)

# Parse options ----

option_list <- list(
    make_option(
        c("--tenxdir", "-t"), action="store",
        type="character",
        dest="tenxdir",
        help="Path to the input 10x matrix directory."),
    make_option(
        c("--fdr", "-s"), action="store",
        type="character",
        dest="fdr",
        default=0.01,
        help="False discovery rate to call cells from ambient noise."),
    make_option(
        c("--outfile", "-s"), action="store",
        type="character",
        dest="outfile",
        help="Description of input option")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")
print(opt)

# Process data ----

# message/done
sce <- read10xCounts(opt$tenxdir, col.names=TRUE)

# message/done
out <- emptyDrops(assay(sce, "counts"))

# write `out` table for the record and diagnostic plots

is.cell <- out$FDR <= 0.01
cat("Number of cells: ")
cat(sum(is.cell, na.rm=TRUE), "\n")

# Check if p-values are lower-bounded by 'npts'
# (increase 'npts' if any Limited==TRUE and Sig==FALSE)
cat("Check if p-values are lower-bounded by 'niters'\n")
print(table(Sig=is.cell, Limited=out$Limited))

# message/done
csvdbTable <- data.frame(
    sample=gsub("([[:alnum:]]+)-count\\/.*", "\\1", opt$tenxdir),
    barcode=colnames(sce),
    emptyDrops=is.cell
)

# message/done
write.table(csvdbTable, opt$outfile, col.names = TRUE, row.names = FALSE, sep = "\t")

# Conclusion ---

message("Completed")
timestamp()
