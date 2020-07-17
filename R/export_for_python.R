## Export data to be used as input for analysis with python
## packages such as scanpy

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
    make_option(c("--counts"), default=FALSE,
                help="Export the raw counts (counts)"),
    make_option(c("--data"), default=FALSE,
                help="Export the normalised data (data)"),
    make_option(c("--scaled"), default=FALSE,
                help="Export the scaled data (data.scaled)"),
    make_option(c("--outdir"), default=".",
                help="the file to which the reduced dimensions will be written")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

## Functions
exportData <- function(seurat_object, slot="counts", outdir=NULL) {
    x <- GetAssayData(s, slot=slot)
    write.table(x, gzfile(file.path(outdir,paste("assay", slot, "tsv.gz", sep="."))),
                quote=FALSE, sep="\t", row.names = FALSE, col.names= FALSE)
    }

exportEmbedding <- function(seurat_object, embedding="PCA", outdir=NULL) {
    x <- Embeddings(object = s, reduction = embedding)
    write.table(x, gzfile(file.path(outdir, paste("embedding", embedding, "tsv.gz", sep="."))),
                quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
    }

exportMetaData <- function(seurat_object, outdir=NULL) {
    x <- seurat_object[[]]
    x$barcode <- Cells(seurat_object)
    write.table(x, gzfile(file.path(outdir, "metadata.tsv.gz")),
                quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
    }


# Read RDS seurat object
if (endsWith(opt$seuratobject, ".rds")) {
   message(sprintf("readRDS: %s", opt$seuratobject))
   s <- readRDS(opt$seuratobject)
} else {
   message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
   stopifnot(require(SeuratDisk))
   s <- LoadH5Seurat(opt$seuratobject)
}

message("export_for_python running with default assay: ", DefaultAssay(s))

# Write out the cell and feature names
message("writing out the cell and feature names")
writeLines(Cells(s), gzfile(file.path(opt$outdir,"barcodes.tsv.gz")))
writeLines(rownames(s), gzfile(file.path(opt$outdir,"features.tsv.gz")))

# Write out embeddings (such as e.g. PCA)
message("Writing matrix of reduced dimensions")
exportEmbedding(s, opt$reductiontype, outdir=opt$outdir)

# Write out the metadata
message("Writing out the metadata")
exportMetaData(s, outdir=opt$outdir)


# Write out significant components
if (opt$usesigcomponents == TRUE) {
  message("Writing vector of significant components")
  comps <- getSigPC(s)
  write.table(comps, file = paste0(opt$outdir, "/sig_comps.tsv"),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

# Write out assay data (such as e.g. raw counts)
if(opt$counts) { exportData(s, "counts", outdir=opt$outdir) }
if(opt$data) { exportData(s, "data", outdir=opt$outdir ) }
if(opt$scaled) { exportData(s, "scale.data", outdir=opt$outdir) }

message("export_for_python.R final default assay: ", DefaultAssay(s))

message("Completed")
