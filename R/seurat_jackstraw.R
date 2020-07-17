## Perform the Seurat JackStraw analysis, which can be useful
## for deciding on the number of PCA components to include in the
## downstrem analysis

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="must specify",
                help="A seurat object with precomputed pca analysis"),
    make_option(c("--numreplicate"), type="integer", default=100,
                help="Number of replicates"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

if (endsWith(opt$seuratobject, ".rds")) {
  message(sprintf("readRDS: %s", opt$seuratobject))
  s <- readRDS(opt$seuratobject)
} else {
  message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
  stopifnot(require(SeuratDisk))
  s <- LoadH5Seurat(opt$seuratobject)
}

## In Macosko et al, we implemented a resampling test inspired by the jackStraw procedure.
## We randomly permute a subset of the data (1% by default) and rerun PCA,
## constructing a 'null distribution' of gene scores, and repeat this procedure. We identify
## 'significant' PCs as those who have a strong enrichment of low p-value genes.

nPCs <- min(dim(s@reductions$pca@cell.embeddings)[2],30)

s <- JackStraw(s, num.replicate=opt$numreplicate,
               num.pc = nPCs)

s <- JackStrawPlot(s, PCs=1:nPCs)

gp <- s@reductions$pca@misc$jackstraw.plot

save_ggplots(paste0(opt$outdir,"/pcaJackStraw"),
           gp,
           width=8,
           height=12)
