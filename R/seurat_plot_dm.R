## Compute reduced diffusion map dimensions and save to a table
## along with the cell metadata

# Libraries ----

stopifnot(
    require(Seurat),
    require(destiny),
    require(plot3D),
    require(tenxutils),
    require(dplyr),
    require(Matrix),
    require(reshape2),
    require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--diffusionmap"), default="dm.tsv.gz",
                help="A seurat object after PCA"),
    make_option(c("--clusterassignments"), default="none",
                help="A tsv files containing the cluster identities"),
    make_option(c("--plotdirvar"), default="diffmapDir",
                help="name of the latex var containing the directory with the plots"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="the directory in which plots will be saved")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


message("read in the diffusion map")
dm <- read.table(opt$diffusionmap, header=T, sep="\t", as.is=T)

message("read in the cluster ids")
#$cluster
clusters <- read.table(opt$clusterassignments, sep="\t", header=T, as.is=T)
rownames(clusters) <- clusters$barcode

dm$cluster <- clusters[dm$barcode,"cluster_id"]

message("drawing the plots")

## make a sensible set of 3D plots.
draw3D <- function(m, phi=30, theta=30, clusters, cols)
{
    scatter3D(m$DC1, m$DC2, m$DC3,
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
        n = length(unique(dm$cluster))
        cols = gg_color_hue(n)

        par(mfrow=c(2,2),
            mai = c(0.1, 0.1, 0.1, 0.1))

        draw3D(dm, phi=30, theta=30, dm$clusters, cols)
        draw3D(dm, phi=30, theta=120, dm$clusters, cols)
        draw3D(dm, phi=30, theta=210, dm$clusters, cols)
        draw3D(dm, phi=30, theta=300, dm$clusters, cols)
    }

plotfilename="diffusion.map"

save_plots(file.path(opt$outdir, plotfilename),
           draw3Dplots,
           width=10,
           height=10)

texCaption <- paste("Diffusion map plots (first 3 dimensions, different rotations) colored by cluster")

tex <- paste(getSubsubsectionTex(texCaption),
             getFigureTex(plotfilename,texCaption,plot_dir_var=opt$plotdirvar),
             sep="\n")

message("Writing out latex snippet")
## write out latex snippet

tex_file <- file.path(opt$outdir,
                      "plot.diffusion.map.tex")

writeTex(tex_file, tex)


message("Completed")
