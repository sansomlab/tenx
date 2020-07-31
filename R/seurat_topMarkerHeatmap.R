## Summarise the marker genes across a set of clusters

# Libraries ----

stopifnot(
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(data.table),
  require(openxlsx),
  require(optparse),
  require(ComplexHeatmap),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--seuratassay"), default="RNA",
                help="The seurat assay to use"),
    make_option(c("--slot"), default="scale.data",
                help="The seurat slot to use for the marker heatmap"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--markers"), default="none",
                help="The table of summarised markers"),
    make_option(c("--subgroup"), default=NULL,
                help="Optional. Name of a column in metadata to add annotation for in the heatmap"),
    make_option(c("--pdf"), default = FALSE,
                help="Create a pdf version of the top marker heatmap"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

outPrefix <- file.path(opt$outdir,"markers.summary")
if(!is.null(opt$subgroup)) { opt$subgroup <- strsplit(opt$subgroup,",")[[1]]}
cat("Running with options:\n")
print(opt)

s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
Idents(s) <- cluster_ids

message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("seurat_summariseMarkers.R running with default assay: ", DefaultAssay(s))

message("Making a heatmap of the top marker genes from each cluster")

markers <- read.table(opt$markers,
                      sep="\t", header=T, as.is=T)

filtered_markers <- data.table(markers[markers$p.adj < 0.1,])

## make a heatmap of the top DE genes.
filtered_markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20

if(!is.null(opt$subgroup))
{
    if(!opt$subgroup %in% colnames(s@meta.data))
    {
        opt$subgroup <- NULL
    }
}


if(max(dim(GetAssayData(s, slot="scale.data"))) > 0 & opt$slot=="scale.data") {
    hm_slot = "scale.data" } else { hm_slot = "data" }

message("using slot: ", hm_slot, " for the top marker heatmap")

mch <- markerComplexHeatmap(s,
                            marker_table=filtered_markers,
                            n_markers=20,
                            slot=hm_slot,
                            cells_use=NULL,
                            row_names_gp=11,
                            sub_group=opt$subgroup)

drawHeatmap <- function()
{
    draw(mch)
}

save_plots(paste(outPrefix,"heatmap", sep="."),
           plot_fn=drawHeatmap,
           width = 7,
           to_pdf = opt$pdf,
           height = 9)


message("seurat_summariseMarkers.R final default assay: ", DefaultAssay(s))

message("completed")
