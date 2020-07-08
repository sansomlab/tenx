# tenx

A collection of Python3 pipelines that call R and Python3 scripts for the analysis of data generated with the 10x Genomics platform. The pipelines are based on 10x's [Cell Ranger pipeline](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) and [DropEst](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1449-6) for mapping and quantitation. 

Downstream analysis currently relies on both the [R Seurat library](https://satijalab.org/seurat/) and [Python Scanpy package](https://scanpy.readthedocs.io/en/stable/), and makes use of many excellent tools from the community including [Scran](https://www.rdocumentation.org/packages/scran/versions/1.0.3), [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html), [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html), [Clustree](https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html), [Destiny](https://bioconductor.org/packages/release/bioc/html/destiny.html) (for diffusion maps), [PHATE](https://www.krishnaswamylab.org/projects/phate), [PAGA](https://github.com/theislab/paga) and [Scvelo](https://scvelo.readthedocs.io/)) for downstream analysis. Automatic export of [UCSC cell browser](https://cells.ucsc.edu/)  instances is also supported.

For geneset over representation analysis the pipelines use a bespoke R package called [gsfisher](https://github.com/sansomlab/gsfisher), which can also be used [interactively to analyse single-cell data](https://github.com/sansomlab/gsfisher/blob/master/vignettes/single_cell_over_representation_analysis.pdf).

The pipelines are in active development, and should be considered "beta" software - please use at your own risk!

# Examples

###  Interferon beta stimulated PBMCs

This example shows how the [Seurat stimulated and control vignette](https://satijalab.org/seurat/v3.1/immune_alignment.html) can be reproduced by the pipeline.

* [Summary Report](https://dl.dropbox.com/s/84x0m9sjdsah8b3/summaryReport.pdf)

* [Aligned example vignette](docs/AlignedExample.md).

### Pancreatic embryogenesis 

This is the [scvelo Bastidas-Ponce et al. dataset](https://scvelo.readthedocs.io/scvelo.datasets.pancreas.html)

* [Summary Report](https://dl.dropbox.com/s/n355kakx6d2jbqp/summaryReport.pdf)
* Full details to follow.

### Microwell-seq Mouse Atlas (240k cells)

Here pipeline_seurat.py was run usin the seurat object provided by the Seurat authors in their [Guided Clustering of the Microwell-seq Mouse Cell Atlas vignette](https://satijalab.org/seurat/v3.1/mca.html). 

* [Summary Report](https://dl.dropbox.com/s/0nyfg5xlsx6u3v1/summaryReport.pdf)
* Full details to follow.


# installation and dependencies

1. [Installation](docs/INSTALL.md)
2. [Dependencies](docs/DEPENDENCIES.md)


# typical workflow

1. Perform mapping, quantification, aggregation and down-sampling using: **`pipelines/pipeline_cellranger.py`**
   * Can be run either from `cellranger mkfastq` or `cellranger aggr` outputs.
   * Samples are mapped and quantitated with `cellranger count`.
   * Aggregation of sample matrices is performed with with `cellranger aggr`.
   * Cells with barcodes shared between cells can be removed (within sequencing batch) to mitigate index hopping.
   * Random down-sampling of the UMI-count matrix is supported.
   * Arbitrary subsets of the aggregated dataset can be generated.

2. Perform downstream analysis using a Seurat based workflow: **`pipelines/pipeline_seurat.py`**
   * This can be run either from count matrices (e.g. `pipeline_cellranger.py` output) or from saved Seurat object(s).
   * Analysis of multiple samples with different parameter combinations can be executed in parallel.
   * Supports testing for genes differently expressed between conditions.
   * Supports finding conserved markers (both between cluster and condition).
   * Support for basic geneset over-enrichment analysis (including of arbitrary "gmt" genesets e.g. from [MSigDB](https://software.broadinstitute.org/gsea/msigdb/)) using [gsfisher](https://github.com/sansomlab/gsfisher).
   * Support for visualising expression of arbitrary lists of genes on violin and UMAP plots.
   * The pipeline includes Clustree, PAGA, ScVelo, Diffusion maps, PHATE maps and SingleR.
   * The pipeline can automatically generate UCSC cell browser instances.
   





