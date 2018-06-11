# tenx

A collection of python3 pipelines and Rscripts for analysing data generated with the 10x Genomics platform. The pipelines are based on 10x's [Cell Ranger pipeline](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) for mapping and quantitation and the [R Seurat package](https://satijalab.org/seurat/) for downstream analysis.

The piplines are in active development, this is an "alpha" release - use at your own risk!


# installation and dependencies

1. [Installation](docs/INSTALL.md)
2. [Dependencies](docs/DEPENDENCIES.md)


# typical workflow

1. Perform mapping, quantitation, aggregation and downsampling using: **`pipelines/pipeline_cellranger.py`**
   * Can be run either from `cellranger mkfastq` or `cellranger aggr` outputs.
   * Samples are mapped and quantiated with `cellranger count`.
   * Aggregation of sample matrices is performed with with `cellranger aggr`.
   * Cells with barcodes shared between cells can be removed (within sequencing batch) to mitigate index hopping.
   * Random downsampling of the umi-count matrix is supported.
   * Arbitrary subsets of the aggregated dataset can be generated.

2. Perform analysis using Seurat: **`pipelines/pipeline_seurat.py`**
   * Can be run either from count matrices (e.g. `pipeline_cellranger.py` output) or from saved Seurat object(s).
   * Analysis of multiple samples with different parameter combinations can be executed in parallel.
   * Supports testing for differences between conditions.
   * Supports finding conserved markers (both between cluster and condition).
   * Support for basic geneset over-enrichment analysis (including of arbitrary "gmt" genesets e.g. from [MSigDB](https://software.broadinstitute.org/gsea/msigdb/)).
   * Support for visualising expression of arbitrary lists of genes on tSNE plots.


# pipeline_seurat.py example.

pipeline_seurat.py was run on the interferon beta stimulated PBMC example dataset from the [Seurat website](https://satijalab.org/seurat/).

The summary report can be download here: [Summary Report](https://dl.dropbox.com/s/67z5xydxvhqdw3p/summaryReport.pdf)

The steps followed, and full output are available here: [aligned example](docs/AlignedExample.md).
