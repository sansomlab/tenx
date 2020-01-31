# pipeline_seurat.py: aligned example

Prior to running the pipeline the [Satija lab vignette](https://satijalab.org/seurat/v3.1/immune_alignment.html) was followed to perform [CCA alignment](https://doi.org/10.1038/nbt.4096) as below. Note that we also scale the data in the RNA assay.

```
# (in R)
# 
library(Seurat) # Version 3.1.1

# devtools::install_github('satijalab/seurat-data')
library(SeuratData)

library(cowplot)

# Work around issue with InstallData("ifnb") 
# (see: https://github.com/satijalab/seurat-data/issues/15)
install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")
library(ifnb.SeuratData)

# load the data
data("ifnb")

ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

# Compute the PCAs on the integrated data
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

# Scale the RNA assay data (needed for visualisations of gene expression levels)
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)

# Switch the default slot back to integrated.
DefaultAssay(immune.combined) <- "integrated"

saveRDS(immune.combined, "immune.combined.rds")
```

pipeline_seurat.py was then configured and run as follows:

```
# (from a bash shell)
# $align_path is the directory containing the "immune.combined.rds" file.
mkdir seurat_aligned_example
cd seurat_aligned_example
mkdir aligned.seurat.dir
ln -s $align_path/immune.aligned.rds aligned.seurat.dir/begin.rds

# generate the configuration file
python $tenx_path/pipelines/pipeline_seurat.py config

# edit the file appropriately, e.g.
emacs -nw pipeline.yml

# run the pipeline
# $tenx_path is the directory containing your clone of the tenx code.
python $tenx_path/pipelines/pipeline_seurat.py make report -v5 -p200
```

The configuration file is available here: [pipeline.yml](https://dl.dropbox.com/s/kvy2r70h9giasie/pipeline.yml)


# Full outputs

The pipeline produces the following files:

1. [Summary Report](https://dl.dropbox.com/s/67z5xydxvhqdw3p/summaryReport.pdf)
2. [Gene expression Report](https://dl.dropbox.com/s/7vq8kxh7kggv7l3/geneExpressionReport.pdf)
3. [Cluster marker genes](https://dl.dropbox.com/s/w0qerus5m2ip7xl/markers.summary.table.xlsx)
4. [Genesets over-enriched in cluster markers](https://dl.dropbox.com/s/l4a2mejov9vfpkr/geneset.analysis.xlsx)
5. [Genes differentially expressed between condition](https://dl.dropbox.com/s/qry6u27l1rxuorx/markers.between.stim.summary.table.xlsx)
6. [Genesets over-enriched in genes differentially expressed between condition](https://dl.dropbox.com/s/nfpunhjgoi0gm4o/geneset.analysis.between.xlsx)
