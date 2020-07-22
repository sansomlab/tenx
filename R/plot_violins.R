## Visualise the expression of given groups of genes
## using horizontal violin plots.

# Libraries ----

stopifnot(
  require(optparse),
  require(ggplot2),
  require(reshape2),
  require(Seurat),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
      c("--genetable"),
      default="none",
      help="A tab-delimited text file containing gene_id (or gene if s@misc$gene_id is not set) and gene_name"
      ),
    make_option(
      c("--seuratobject"),
      default="none",
      help="The seurat object (e.g. begin.rds)"
    ),
    make_option(
      c("--seuratassay"),
      default="RNA",
      help="The assay to set as default before GetAssayData is called"
    ),
    make_option(
      c("--clusterids"),
      default="none",
      help="The rdsfile containing the clusterids"
    ),
    make_option(
      c("--plotdirvar"),
      default="clusterMarkerTSNEPlotsDir",
      help="The name of the latex var specifying the location of the plots"
      ),
    make_option(
      c("--outprefix"),
      default="seurat.out.dir",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


## read in the table containing the genes to visualise
genes <- read.table(opt$genetable,
                    header=TRUE,
                    sep="\t",
                    as.is=TRUE)


## read in the raw count matrix
s <- loadSeurat(path=opt$seuratobject)

Idents(s) <- readRDS(opt$clusterids)

## set the default assay
message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay


message("plot_violins.R running with default assay: ", DefaultAssay(s))

if("gene_id" %in% colnames(s@misc))
{
    ## map to the seurat ids via the ensembl_ids...
    id_col = "gene_id"
    map <- data.frame(s@misc[s@misc$gene_id %in% genes$gene_id,])

    #subset the map to genes with expression data.
    map <- map[map$seurat_id %in% rownames(x = s),]
    rownames(map) <- map$gene_id
    genes <- genes[genes$gene_id %in% rownames(map),]
    genes$seurat_id <- as.vector(map[genes$gene_id,"seurat_id"])

    message("mapped to seurat gene via ensembl_id")

} else {
    ## map directly using gene symbols
    id_col = "gene"
    genes <- genes[genes$gene %in% rownames(x = s),]
    genes$seurat_id <- genes$gene

    message("mapped gene symbol directly to seurat gene...")
}


tex = ""

for(group in unique(genes$group))
{
    tmp <- genes[genes$group==group,]


    if("plot_name" %in% colnames(genes))
    {
        plot_names <- make.unique(paste(tmp$plot_name," (",tmp$gene,")",sep=""))
    } else {
        plot_names <- make.unique(tmp$gene)
    }
    tmp$plot_name <- plot_names

    rownames(tmp) <- tmp$seurat_id

    message("prepared the gene table for group: ", group)


    print(rownames(tmp))
    print(rownames(tmp)[rownames(tmp) %in% rownames(x = s)])

    data <- as.data.frame(as.matrix(GetAssayData(object = s, slot="data")[rownames(tmp),]))
    rownames(data) <- tmp$plot_name

    data$gene <- as.vector(rownames(data))

    ggData <- melt(data,id.vars="gene")

    ggData$cluster <- as.numeric(as.vector(Idents(s)[ggData$variable]))

    message("prepared the expression data")

    cluster_uids <- unique(Idents(s))

    colors <- gg_color_hue(length(cluster_uids))
    names(colors) <- cluster_uids[order(as.numeric(cluster_uids))]

    ncol <- 10

    ggGrob <- makeViolins(ggData,
                      title=group,
                      ncol=ncol,
                      xlab="normalised gene expression",
                      clusters=NULL,
                      group=NULL,
                      colors=colors,
                      alpha=1)

    plot_path = paste(opt$outprefix, group, sep=".")

    save_ggplots(plot_path,
                ggGrob,
                width=8,
                height=max(2.5, (ceiling(nrow(tmp)/ncol)*2.5))
                )

    message("saved the plots for group ", group)

    texCaption <- paste(group," genes")
    tex <- c(tex, getFigureTex(basename(plot_path),
                               texCaption,
                               plot_dir_var=opt$plotdirvar))
}

message("Completed plotting.")

tex_file <- file.path(paste(opt$outprefix,"tex",
                            sep="."))

writeTex(tex_file, tex)

message("plot_violins.R final default assay: ", DefaultAssay(s))

message("completed")
