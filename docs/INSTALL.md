# Installation

Check the list of [dependencies](DEPENDENCIES.md)

1. Install the cgat-core pipeline system following the instructions here: https://github.com/cgat-developers/cgat-core/
2. Install the [gsfisher R library](https://github.com/sansomlab/gsfisher), e.g.
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
5. Install velocyto using e.g.
```
   R> library(devtools)
   R> install_github("velocyto-team/velocyto.R")
```
6. Install umap using
```
   $> pip install umap-learn
```