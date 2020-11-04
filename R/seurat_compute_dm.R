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
    make_option(c("--reductiontype"), default="pca",
                help="Name of dimensional reduction technique to use for building the diffusion map. (e.g. 'pca', 'ica')"),
    make_option(c("--outfile"), default="dm.tsv",
                help="the file to which the diffusion map coordinates will be written")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


message("readRDS")
s <- readRDS(opt$seuratobject)

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

dm_coord$barcode <- row.names(dm_coord)



## save the annotated coordinates
write.table(dm_coord, gzfile(opt$outfile),
            sep="\t", quote=FALSE, row.names=FALSE)




message("seurat_dm.R final default assay: ", DefaultAssay(s))

message("Completed")
