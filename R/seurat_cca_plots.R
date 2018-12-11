## Make the CCA plots

# Libraries ----

stopifnot(
    require(Seurat),
    require(tenxutils),
    require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.rds",
                help="A single seurat object containing all the samples to be aligned"),
    make_option(c("--metavar"), default=NULL,
                help="The meta data variable that stratifies the samples to be aligned"),
    make_option(c("--strata"), default="all",
                help="the levels of the meta var that should be included in the alignment"),
    make_option(c("--numccs"), type="integer", default=1,
                help="number of ccs to evaluate for the bicor plot."),
    make_option(c("--outprefix"), default="cca.plots",
                help="the name for the output file")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


cca <- readRDS(opt$seuratobject)


message("doing the biocor plot")
## Make the metagene bicorrelation plot
gp <- MetageneBicorPlot(cca,
                        grouping.var = opt$metavar,
                        dims.eval = 1:opt$numccs,
                        display.progress = FALSE)

save_ggplots(paste0(opt$outprefix,".bicor"),
             gp=gp,
             width=6,
             height=4)


message("drawing the component heatmaps")
## Visualise the components

drawHeatmap <- function()
    {
        DimHeatmap(object = cca,
                 reduction.type = "cca",
                 cells.use = 500,
                 dim.use = 1:12,
                 do.balanced = TRUE)
    }

save_plots(paste0(opt$outprefix,".heatmaps"),
             plot_fn=drawHeatmap,
             width=5,
             height=8)

message("completed")
