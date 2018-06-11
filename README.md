# tenx

A collection of pipelines and Rscripts for analysing data generated with the 10x Genomics platform. The pipelines are based on 10x's [Cell Ranger pipeline](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) for mapping and quantitation and the [R Seurat package](https://satijalab.org/seurat/) for downstream analysis.

The piplines are in active development, this is an "alpha" release - use at your own risk!


# installation and dependencies

1. [Installation](docs/INSTALL.md)
2. [Dependencies](docs/DEPENDENCIES.md)


# typical workflow

1. Run the 10x Cell Ranger pipeline and perform downsampling using: **pipelines/pipeline_cellranger.py**
2. Perform analysis using Seurat: **pipelines/pipeline_seurat.py**


# pipeline_seurat.py example.

pipeline_seurat.py was run on the interferon beta stimulated PBMC example dataset from the [Seurat websit](https://satijalab.org/seurat/). Prior to running the pipeline the [vignette](https://satijalab.org/seurat/immune_alignment.html) was followed to perform [CCA alignmet](https://doi.org/10.1038/nbt.4096).

The summary report can be download here: [Summary Report](https://dl.dropbox.com/s/67z5xydxvhqdw3p/summaryReport.pdf)

The steps followed, and full output are available here: [aligned example](docs/AlignedExample.md).
