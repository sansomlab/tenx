# Script to generate input files for cellbrowser

# Libraries ----
stopifnot(
  require(optparse),
  require(Seurat),
  require(gridExtra),
  require(openxlsx),
  require(tidyr),
  require(dplyr),
  require(tidyselect),
  require(reshape2),
  require(tenxutils)
)

# Options ----

option_list <- list(
  make_option(
    c("--seurat_path"),
    help="Location of the begin.rds (has to include annotation in misc slot)."
  ),
  make_option(
    c("--runspecs"),
    help="Clustering runspecs as folder name string"
  ),
  make_option(
    c("--outdir"),
    default="seurat.out.dir",
    help="Location for outputs files. Must exist."
  ),
  make_option(
    c("--add_fdg"),
    default=NULL,
    help="Whether to add FDG coordinates to the cellbrowser"
  ),
  make_option(
    c("--add_phate"),
    default=NULL,
    help="Whether to add PHATE coordinates to the cellbrowser"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


cat("Read Seurat object ... \n")
seurat_obj = readRDS(file.path(opt$seurat_path, "begin.rds"))

# read in annotation to add gene names
cat("Extract gene annotations ... \n")
annotation <- seurat_obj@misc


## Load UMAP coordinates
cat("Load UMAP coordinates ... \n")
runspecs = data.frame(all = opt$runspecs) %>% separate(all, c("comp", "res", "mindist"), sep = "_")
data_selected = read.table(file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                                     "umap.dir", paste0("umap.", runspecs$mindist,".tsv.gz")),
                           header=TRUE, as.is=TRUE, sep="\t")
output = data_selected[,c("barcode", "UMAP_1", "UMAP_2")]
colnames(output) = c("cellId", "x", "y")
write.table(output,file.path(opt$outdir, "UMAP.tsv"), quote = FALSE, row.names = FALSE,
            sep = "\t")

if (!is.null(opt$add_fdg)){
  ## Load Force-directed graph coordinates
  cat("Load Force-directed graph coordinates ... \n")
  infile_gzip = gzfile(file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                                 paste0("cluster.", runspecs$res, ".dir"),  "paga.dir", "paga_init_fa2.tsv.gz"))
  output = read.table(infile_gzip, header=TRUE, as.is=TRUE, sep="\t")
  output = output[,c("barcode", "FA1", "FA2")]
  colnames(output) = c("cellId", "x", "y")
  write.table(output,file.path(opt$outdir, "FA.tsv"), quote = FALSE, row.names = FALSE,
              sep = "\t")
}

if (!is.null(opt$add_phate)){
  cat("Load PHATE coordinates ... \n")
  infile_gzip = gzfile(file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                                 "phate.dir", "phate.tsv.gz"))
  output = read.table(infile_gzip, header=TRUE, as.is=TRUE, sep="\t")
  output = output[,c("barcode", "PHATE1", "PHATE2")]
  colnames(output) = c("cellId", "x", "y")
  write.table(output,file.path(opt$outdir, "PHATE.tsv"), quote = FALSE, row.names = FALSE,
              sep = "\t")
}


## get metadata for cluster name
cat("Process cluster ids for cells ... \n")
infile_cluster_assignments = gzfile(file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                                               paste0("cluster.", runspecs$res, ".dir"), "cluster_assignments.tsv.gz"))
output = read.table(infile_cluster_assignments, header=TRUE, as.is=TRUE, sep="\t")
colnames(output) = c("barcode", "cluster")
print(table(output$cluster))

cat("Process other metadata ... \n")
meta = read.table(gzfile(file.path(opt$seurat_path, "metadata.tsv.gz")), 
                           header=TRUE, as.is=TRUE, sep="\t")
output = merge(output, meta, on='barcode')
write.table(output,file.path(opt$outdir, "meta.tsv"), quote = FALSE, row.names = FALSE,
            sep = "\t")

## write out.tsv of expression for all included cells
cat("Prepare expression matrix ... \n")
expr_data = as.data.frame(as.matrix(GetAssayData(seurat_obj, slot = "data")))
expr_data = cbind(data.frame(gene=rownames(expr_data)),expr_data)
print(dim(expr_data))
out = gzfile(file.path(opt$outdir, "exprMatrix.tsv.gz"), "wt")
write.table(expr_data, out, row.names = FALSE, quote = FALSE,
            sep = "\t", col.names = TRUE)
close(out)
cat("Finished writing expression matrix ... \n")

## set the cell colors according to cluster
cat("Read in the colors to match colors in other visualisations ... \n")

cols = read.csv(file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                          paste0("cluster.", runspecs$res, ".dir"), "cluster_colors.tsv"),
                  header = FALSE, stringsAsFactors = FALSE)
clusters = sort(unique(as.numeric(output$cluster)))
out_color = data.frame(name = clusters,
                       color = cols,
                       stringsAsFactors = FALSE)
out_color$color = substr(out_color$V1, 2, nchar(out_color$V1)[1])
out_color$V1 = NULL
write.table(out_color, file.path(opt$outdir,"colors.tsv"), quote = FALSE,
            sep = "\t", row.names = FALSE)

## add marker genes
cat("Prepare marker genes to be shown ... \n")
print(file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                paste0("cluster.", runspecs$res, ".dir"), "cluster.markers.dir","markers.summary.table.xlsx"))
data_selected = read.xlsx(xlsxFile = file.path(opt$seurat_path, paste0("components.", runspecs$comp, ".dir"), 
                                               paste0("cluster.", runspecs$res, ".dir"),
                                              "cluster.markers.dir","markers.summary.table.xlsx"))
output = data_selected[,c("gene","gene_id", "cluster","avg_logFC","p.adj")]
output = output %>% group_by(cluster) %>% dplyr::arrange(desc(avg_logFC)) %>% do(head(.,n=20)) %>% ungroup()
output$celltype = "top20_marker_Seurat_cluster"
output$cluster_marker = paste(output$celltype, output$cluster, sep="_")
output = output %>% dplyr::select(-celltype)
colnames(output) = c("gene","gene_id","cluster","avg_logFC","p_adjusted","celltype_marker")
output = output[,c("cluster","gene","p_adjusted","avg_logFC","celltype_marker")]
write.table(output, file.path(opt$outdir, "markers.tsv"), sep = "\t",
          quote = FALSE, row.names = FALSE)

cat("Completed")