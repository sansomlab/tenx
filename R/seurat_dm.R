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
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--annotation"), default="none",
                help="A file containing the mapping of gene_id, gene_name and seurat_id"),
    make_option(c("--usegenes"), default=FALSE,
                help="use variable genes instead of reduced dimensions"),
    make_option(c("--usesigcomponents"), default=FALSE,
                help="use significant principle component"),
    make_option(c("--components"), type="integer", default=10,
                help="number of principle components to use"),
    make_option(c("--maxdim"), type="integer", default=10,
                help="the choice of metric used to measure distance in the input space"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use for building the diffusion map. (e.g. 'pca', 'ica')"),
    make_option(c("--plotdirvar"), default="diffmapDir",
                help="name of the latex var containing the directory with the plots"),
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
Idents(s) <- cluster_ids

message("seurat_dm.R running with default assay: ", DefaultAssay(s))

## run the diffusion map algorithm
if(opt$usegenes)
{
message("Running diffusion with the variable genes")

diff.data <- t(GetAssayData(s, slot = "scale.data")[VariableFeatures(s),])

} else {

    ## get the principle components to use
    if(opt$usesigcomponents)
    {
        comps <- getSigPC(s)
        message("using the following pcas:")
        print(comps)
    } else {
        comps <- 1:as.numeric(opt$components)
    }

    message("Running diffussion with reduced dimensions")

    diff.data <- Embeddings(s, reduction=opt$reductiontype)[,comps]
}

message("making the diffusion map")

print(dim(diff.data))

## make the diffusion map
diff_map <- DiffusionMap(diff.data)

message("saving the dm coordinates")

## extract the DM coordinates from the seurat object
dm_coord <- as.data.frame(diff_map@eigenvectors)

print(head(dm_coord))
print(dim(dm_coord))
print(head(Idents(s)))

rownames(dm_coord) <- Cells(s)
dm_coord$cluster <- Idents(s)

dm_coord <- merge(dm_coord, s[[]], by=0)

dm_coord$barcode <- row.names(dm_coord)

## save the annotated coordinates
write.table(dm_coord, gzfile(opt$outfile),
            sep="\t", quote=FALSE, row.names=FALSE)


message("drawing the plots")

## make a sensible set of 3D plots.
draw3D <- function(m, phi=30, theta=30, clusters, cols)
{
    scatter3D(m@eigenvectors[,1],
              m@eigenvectors[,2],
              m@eigenvectors[,3],
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
        n = length(unique(Idents(s)))
        cols = gg_color_hue(n)
        clusters <- as.numeric(as.vector(Idents(s)))

        par(mfrow=c(2,2),
            mai = c(0.1, 0.1, 0.1, 0.1))

        draw3D(diff_map, phi=30, theta=30, clusters, cols)
        draw3D(diff_map, phi=30, theta=120, clusters, cols)
        draw3D(diff_map, phi=30, theta=210, clusters, cols)
        draw3D(diff_map, phi=30, theta=300, clusters, cols)
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


message("seurat_dm.R final default assay: ", DefaultAssay(s))

message("Completed")
