# pipeline_seurat.py: aligned example

```
# (in R)
# Follow step from: https://satijalab.org/seurat/immune_alignment.html
# up to and including
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:20)

# save the combined object for pipeline_seurat
saveRDS(immune.combined, "immune.combined.rds")
```

The pipeline was then run as follows:
```
# (from a bash shell)
mkdir seurat_aligned_example
cd seurat_aligned_example
mkdir aligned.seurat.dir
ln -s $align_path/immune.aligned.rds aligned.seurat.dir/begin.rds

# generate the configuration file
python $tenx_path/pipelines/pipeline_seurat.py config

# edit the file appropriately, e.g.
emacs -nw pipeline_seurat.py

## run the pipeline
python $tenx_path/pipelines/pipeline_seurat.py -v5 -p200 make report
```

# Expected outputs

[Summary Report](https://dl.dropbox.com/s/67z5xydxvhqdw3p/summaryReport.pdf)
[Gene expression Report](https://dl.dropbox.com/s/7vq8kxh7kggv7l3/geneExpressionReport.pdf)
[Cluster marker genes](https://dl.dropbox.com/s/w0qerus5m2ip7xl/markers.summary.table.xlsx)
[Genesets over-enriched in cluster markers](https://dl.dropbox.com/s/l4a2mejov9vfpkr/geneset.analysis.xlsx)
[Genes differentially expressed between condition](https://dl.dropbox.com/s/qry6u27l1rxuorx/markers.between.stim.summary.table.xlsx)
[Genesets over-enriched in genes differentially expressed between condition](https://dl.dropbox.com/s/nfpunhjgoi0gm4o/geneset.analysis.between.xlsx)
[Configuration file](https://dl.dropbox.com/s/kvy2r70h9giasie/pipeline.yml)