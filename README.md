# tenx

A collection of pipelines and Rscripts for analysing data generated with the 10x Genomics platform. The pipelines are based on 10x's [Cell Ranger pipeline](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) for mapping and quantitation and the [R Seurat package](https://satijalab.org/seurat/) for downstream analysis.

# Dependencies

* Python3
* The [cgat-core pipeline system](https://github.com/cgat-developers/cgat-core/) and a compatible queue manager
* R >= 3.42
* R libraries:
  * Seurat
  * cellrangerRkit
  * ComplexHeatmap
  * data.table
  * DESeq2
  * dplyr
  * ggplot2
  * ggrepel
  * gplots
  * RColorBrewer
  * R.utils
  * grid
  * gridExtra
  * Matrix
  * openxlsx
  * optparse
  * reshape2
  * S4Vectors
  * xtable
  * gsfisher (see below)
  * tenxutils (see below)
* Latex
  
# Installation

1. Install the cgat-core pipeline sytem following the instructions here: https://github.com/cgat-developers/cgat-core/
2. Install the gsfisher R library, e.g.
```
   R> library(devtools)
   R> install_github("sansomlab/gsfisher")
```
3. Clone the tenx repository, e.g.
```
   $> git clone https://github.com/sansomlab/tenx.git
```
4. Install the tenxutils R library e.g.
```
   $> cd tenx
   $> R CMD INSTALL --no-multiarch --with-keep.source tenxutils
```
# Typical workflow

1. Run the 10x Cell Ranger pipeline and perform downsampling using: **pipelines/pipeline_cellranger.py**
2. Perform analysis using Seurat: **pipelines/pipeline_seurat.py**
