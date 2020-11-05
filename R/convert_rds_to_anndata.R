## Export data to be used as input for analysis with python
## packages such as scanpy

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(SeuratDisk),
  require(tenxutils)
)

# Options ----

option_list <- list(
  make_option(c("--seuratobject"), default="begin.Robj",
              help="A seurat object after dimensionality reduction")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

output_file = gsub(".h5seurat", ".h5ad", opt$seuratobject)

Convert(opt$seuratobject, dest = output_file, overwrite = TRUE, verbose = TRUE) 

print("Completed")
