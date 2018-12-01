## Compute reduced diffusion map dimensions and save to a table
## along with the cell metadata

# Libraries ----

stopifnot(
    require(Seurat),
    require(plot3D),
    require(tenxutils),
    require(dplyr),
    require(Matrix),
    require(reshape2),
    require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--annotation"), default="none",
                help="A file containing the mapping of gene_id, gene_name and seurat_id"),
    make_option(c("--components"), type="integer", default=10,
                help="number of principle components to use"),
    make_option(c("--maxdim"), type="integer", default=20,
                help="the choice of metric used to measure distance in the input space"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use for building the diffusion map. (e.g. 'pca', 'ica')"),
    make_option(c("--outfile"), default="dm.txt",
                help="the file to which the diffusion map coordinates will be written"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="the directory in which plots will be saved")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


message("readRDS")
s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
s@ident <- cluster_ids


## run the diffusion map algorithm
message("RunUMAP")
s <- RunDiffusion(object=s,
             reduction.use = opt$reductiontype,
             max.dim=opt$maxdim)

## extract the DM coordinates from the seurat object
dm <- as.data.frame(s@dr$dm@cell.embeddings)
dm$cluster <- s@ident

plot_data <- merge(dm, s@meta.data, by=0)

rownames(plot_data) <- plot_data$Row.names
plot_data$Row.names <- NULL

plot_data$barcode <- row.names(plot_data)

## save the annotated tSNE data frame
write.table(plot_data, opt$outfile,
            sep="\t", quote=FALSE, row.names=FALSE)


## make a sensible set of 3D plots.

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

draw3D <- function(m, phi=30, theta=30, clusters, cols)
{
    scatter3D(m@cell.embeddings[,1],
              m@cell.embeddings[,2],
              m@cell.embeddings[,3],
              bty = "b2",
              col=cols,
              colkey=F,
              phi = phi,
              theta = theta,
              expand=1,
              pch = 19,
              alpha = 0.4,
              colvar=clusters)
}

draw3Dplots <- function()
    {
        n = length(unique(s@ident))
        cols = gg_color_hue(n)
        clusters <- as.numeric(as.vector(s@ident))

        m <- s@dr$dm

        par(mfrow=c(2,2),
            mai = c(0.1, 0.1, 0.1, 0.1))

        draw3D(m, phi=30, theta=30, clusters, cols)
        draw3D(m, phi=30, theta=120, clusters, cols)
        draw3D(m, phi=30, theta=210, clusters, cols)
        draw3D(m, phi=30, theta=300, clusters, cols)
    }

plotfilename="diffusion.map"

save_plots(file.path(opt$outdir, plotfilename),
           draw3Dplots,
           width=10,
           height=10)

texCaption <- paste("Diffusion map plots (first 3 dimensions, different rotations) colored by cluster")



tex <- paste(getSubsubsectionTex(texCaption),
             getFigureTex(plotfilename,texCaption),
             sep="\n")

print("Writing out latex snippet")
## write out latex snippet

tex_file <- file.path(opt$outdir,
                      "plot.diffusion.map.tex")

writeTex(tex_file, tex)



message("Completed")
