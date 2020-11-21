# Installation

Check the list of [dependencies](DEPENDENCIES.md) for more information.

1. Install the cgat-core pipeline system following the instructions here: https://github.com/cgat-developers/cgat-core/

2. Clone the tenx repository, e.g.
```
   $> git clone https://github.com/sansomlab/tenx.git
```

3. In the same virtual or conda environment as cgat-core install the required python packages:
```
   $> pip install -r tenx/python/requirements.txt
```

4. To install the required R packages (R>4.0.0, bioconductor and devtools are prerequiste):
```
   $> Rscript tenx/R/install.packages.R
```
