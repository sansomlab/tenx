## Export data to be used as input for paga

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after dimensionality reduction"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction slot to write (e.g. 'pca', 'integratedreduced')"),
    make_option(c("--usesigcomponents"), default=TRUE,
                help="Whether or not the pipeline is using significant components"),
    make_option(c("--outdir"), default=".",
                help="the file to which the reduced dimensions will be written")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Read RDS seurat object
message("readRDS")
s <- readRDS(opt$seuratobject)
message("paga_prepare_input.R running with default assay: ", DefaultAssay(s))

# Write matrix of reduced dimensions
message("Writing matrix of reduced dimensions")
rd <- Embeddings(object = s, reduction = opt$reductiontype)
out = gzfile(file.path(opt$outdir,"/reduced_dims.tsv.gz"), "wt")
write.table(rd, out,  quote = FALSE, sep = "\t", col.names = TRUE)
close(out)

# Print out significant com
if (opt$usesigcomponents == TRUE) {
  message("Writing vector of significant components")
  comps <- getSigPC(s)
  write.table(comps, file = paste0(opt$outdir, "/sig_comps.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

message("paga_prepare_input.R final default assay: ", DefaultAssay(s))

message("Completed")
