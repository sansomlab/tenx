# pipeline_seurat.py: aligned example

Prior to running the pipeline the [Satija lab vignette](https://satijalab.org/seurat/immune_alignment.html) was followed to perform [CCA alignment](https://doi.org/10.1038/nbt.4096).

```
# (in R)
# Follow step from: https://satijalab.org/seurat/immune_alignment.html
# up to and including
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)

# save the combined object for pipeline_seurat
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
python $tenx_path/pipelines/pipeline_seurat.py -v5 -p200 make report
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
