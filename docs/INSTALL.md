# Installation

Check the list of [dependencies][DEPENDENCIES.md]

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
