stopifnot(
  require(optparse),
  require(tenxutils),
  require(ggplot2),
  require(reshape2),
  require(dplyr),
  require(colormap),
  require(circlize),
  require(Seurat),
  require(ComplexHeatmap)
)


# Options ----

option_list <- list(
  make_option(
    c("--markers"),
    default="none",
    help="A table containing the marker information"),
  make_option(
    c("--clusterids"),
    default="none",
    help="A an rds object containing the clusters"),
  make_option(
    c("--seuratobject"),
    default="none",
    help="The seurat object (typically begin.rds)"
  ),
  make_option(
    c("--seuratassay"),
    default="RNA",
    help="The seurat assay from which the counts will be retrieved",
  ),
  make_option(
    c("--slot"),
    default="data",
    help="the slot to use for the heatmap"
  ),
  make_option(
    c("--group"),
    default=NULL,
    help="A column in the cell metadata used to define the groups"
  ),
  make_option(
    c("--cluster"),
    default=0,
    help="The cluster to generate the plots for"
  ),
  make_option(
    c("--rdimstable"),
    default=NULL,
    help="the reduced dimension coordinates"
  ),
  make_option(
    c("--rdim1"),
    default="UMAP_1",
    help="Rdims x coordinate column"
  ),
  make_option(
    c("--rdim2"),
    default="UMAP_2",
    help='Rdims y coordinate'
  ),
  make_option(
    c("--ngenes"),
    default=24,
    type="integer",
    help='number of genes to pull, currently ignored...'
  ),

  make_option(
    c("--outdir"),
    default="seurat.out.dir",
    help="outdir"),
  make_option(
    c("--pdf"),
    default=FALSE,
    help="outdir")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

outprefix = file.path(opt$outdir,paste("cluster",opt$cluster,sep="."))

# rundir <- "/gfs/work/ssansom/combat/fulldepth/29072020/f_others/"
# mf <- "harmony.seurat.dir/components.50.dir/cluster.1.dir/cluster.markers.dir/markers.summary.table.tsv.gz"
# s <- readRDS(file.path(rundir,"harmony.seurat.dir/begin.rds"))
# cluster_ids <- readRDS(file.path(rundir,"harmony.seurat.dir/components.50.dir/cluster.1.dir/cluster_ids.rds"))
# Idents(s) <- cluster_ids

# read in the Seurat object
s <- readRDS(opt$seuratobject)
Idents(s) <- readRDS(opt$clusterids)

message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay


# read in the markers
x <- read.table(opt$markers, sep="\t", header=T, as.is=T)

# read in the rdims table
rdims <- read.table(opt$rdimstable,
                   sep="\t", header=T,  as.is = T)

# Filter to significant markers for the cluster of interest
x <- x[x$cluster==as.numeric(opt$cluster) & x$p.adj<0.1,]

n_select = 16
# pull out by average and minimum log fold change
x %>% top_n(n_select, avg_log2FC) -> top_by_avg_log2FC

n_select = 8
x[!x$gene %in% top_by_avg_log22FC$gene,] %>%
  top_n(n_select, min_logFC) -> top_by_min_logFC

x <- x[x$gene %in% c(top_by_avg_log2FC$gene,
                     top_by_min_logFC$gene),]

# marker gene heatmap.

if(nrow(x) > 0) {

message("making the heatmap")
mch <- markerComplexHeatmap(s,
                            marker_table=x,
                            n_markers=24,
                            cells_use=NULL,
                            row_names_gp=8,
                            slot="data",
                        #   priority="min_logFC",
                            sub_group=opt$group)

drawHeatmap <- function() { draw(mch) }

save_plots(paste(outprefix,"heatmap", sep="."),
           plot_fn=drawHeatmap,
           width = 7,
           to_pdf = FALSE, #opt$pdf,
           height = 2)

message("making the expression dotplots")
# expression plots
gp <- expressionPlots(s,
                x$gene,
                rdims, x=opt$rdim1, y=opt$rdim2,
                ncol = 4, pch = ".", point_size = 1,
                max_quantile = 0.9)

gp <- gp + theme(panel.spacing.x=unit(0, "lines") ,
                 panel.spacing.y=unit(0, "lines"))

gp <- gp +theme(strip.text = element_text(size = 6, margin = margin()))

## save the plots
save_ggplots(paste(outprefix,"rdims",sep="."),
             gp,
             width=7,
             height=9,
             to_pdf=opt$pdf)

message("making the violin plots")
# make the violin plots
gg_grob <- plotHorizontalViolins(s,
                       x$gene,
                       clusters=NULL,
                       title=NULL,
                       ncol=12,
                       xlab="normalised expression level",
                       group=NULL,
                       colors=NULL,
                       alpha=1,
                       plot=FALSE)


do_plot <- function() { plot(gg_grob)}

do_plot()

nclusters <- length(unique(Idents(s)))

height = min(nclusters/20 * 5,10)

## save the plots
save_plots(paste(outprefix,"violins",sep="."),
             plot_fn=do_plot,
             width=7,
             height=height,
             to_pdf=opt$pdf)


message("making the split dotplot")
# dot plots.
s <- subset(s, cells=names(Idents(s)[Idents(s)==opt$cluster]))

# sort out colors.

## sort out the colors
if(!is.null(opt$group))
{
    fvar <- opt$group
    ngroups <- length(unique(s[[fvar]][[fvar]]))
} else {
    fvar <- NULL
    ngroups <- 1
}

cm_palette <- colormap(colormap = colormaps$portland,
                       nshade = ngroups, alpha=0.6)

message("building the plot")

gp <- gp + scale_fill_manual(values=cm_palette)

gp <- DotPlot(s,
          features=x$gene,
          cols=cm_palette,
          split.by=fvar)

gp <- gp + ylab(fvar) + xlab("gene")
gp <- gp + theme_bw(base_size=8)
gp <- gp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

message("saving the plot")
## save the plots
save_ggplots(paste(outprefix,"dotplot",sep="."),
             gp,
             width=7,
             height=3,
             to_pdf=opt$pdf)


}
message("seurat_summariseMarkers.R final default assay: ", DefaultAssay(s))
message("completed")
