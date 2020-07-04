## Title ----
##
## Visualisation of factors on reduced dimension (2D) plots
##
## Description ----
##
## Visualise the single-cells by their coordinates on reduced dimensions
## Cells are colored by cluster, and by the barcode fields
## which typically represent the sequence batch, sample
## and aggregation id.
##
## Details ----
##
##
## Usage ----
##
## See options.

# Libraries ----

stopifnot(
  require(ggplot2),
  require(reshape2),
  require(optparse),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
      c("--table"),
      default="none",
      help="A table containing the reduced coordinates and phenotype information"
    ),
    make_option(
        c("--metadata"),
        default="none",
        help="A table containing phenotype information"),
    make_option(
      c("--method"),
      default="tSNE",
      help="Normally the type of dimension reduction"
      ),
    make_option(
      c("--rdim1"),
      default="tSNE_1",
      help="The name of the column for reduced dimension one"
      ),
    make_option(
      c("--rdim2"),
      default="tSNE_2",
      help="The name of the column for reduced dimension two"
      ),
    make_option(
      c("--shapefactor"),
      default="none",
      help="A column in the cell metadata to use for deterimining the shape of points on the tSNE"
      ),
    make_option(
      c("--colorfactors"),
      default="none",
      help="Column(s) in the cell metadata to use for deterimining the color of points on the tSNE. One plot will be made per color factor."
    ),
    make_option(
      c("--plotdirvar"),
      default="tsneDir",
      help="latex var holding location of plots"
      ),
    make_option(
      c("--pointsize"),
      default=0.5,
      help="The point size for the tSNE plots"
    ),
    make_option(
      c("--pointpch"),
      default="16",
      help="The point pch for the tSNE plots (if no shapefactor)"
      ),
    make_option(
      c("--pointalpha"),
      default=0.8,
      help="The alpha setting for the points on the tSNE plots"
    ),
    make_option(
      c("--pdf"),
      default=FALSE,
      help="Produce pdf versions of the plots"
      ),
    make_option(
      c("--outdir"),
      default="seurat.out.dir",
      help="outdir"
      )
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


##
message("Reading in rdims table")
plot_data <- read.table(opt$table, sep="\t", header=TRUE)
rownames(plot_data) <- plot_data$barcode

if(opt$metadata!="none")
{
    meta_data <- read.table(opt$metadata, sep="\t", header=TRUE)
    rownames(meta_data) <- meta_data$barcode
    meta_data$barcode <- NULL
    meta_cols <- colnames(meta_data)

    if(length(intersect(rownames(plot_data), rownames(meta_data))) < length(rownames(meta_data)))
    {
        stop("Not all cell barcodes are present in the given metadata table")
    } else {

        meta_data <- meta_data[rownames(plot_data),]
        plot_cols <- colnames(plot_data)
        plot_data <- cbind(plot_data, meta_data)
        colnames(plot_data) <- c(plot_cols, meta_cols)
    }
}

if ("cluster" %in% colnames(plot_data)){
    cluster_col = "cluster"
} else if ("cluster_id" %in% colnames(plot_data))
{
    cluster_col = "cluster_id"
} else
{
    cluster_col = NULL
}


if(!is.null(cluster_col))
{
  plot_data[[cluster_col]] <- factor(plot_data[[cluster_col]], levels=sort(as.numeric(unique(plot_data[[cluster_col]]))))
}

color_vars <- strsplit(opt$colorfactors,",")[[1]]
print(color_vars)
tex = ""

print("Making tSNE plots colored by each of the factor variables")
## Make one whole-page tSNE plot per color variable
for(color_var in color_vars)
{
    print(paste("Making",color_var,"tSNE plot"))


    if(color_var==cluster_col)
    {
        clust_levels = unique(plot_data[[cluster_col]])

        print(head(plot_data))
        message("computing cluster centers")
        centers = data.frame(row.names=clust_levels, "cluster"=clust_levels,
                             "x"=1,
                             "y"=1)

        print(head(centers))
        for(clust in clust_levels)
        {
            centers[clust,"x"] = median(plot_data[plot_data[[cluster_col]]==clust,opt$rdim1])
            centers[clust,"y"] = median(plot_data[plot_data[[cluster_col]]==clust,opt$rdim2])
        }
        message("computed cluster centers")
        print(head(centers))

    }


    ## If a variable comprises only integers, coerce it to a character vector

    numeric = FALSE
    if(is.numeric(plot_data[[color_var]]))
    {
        if((all(plot_data[[color_var]] == round(plot_data[[color_var]]))
           & length(unique(plot_data[[color_var]])) < 50) || color_var == cluster_col)
        {
            plot_data[[color_var]] <- as.character(plot_data[[color_var]])
        }
        else
        {
            numeric = TRUE
        }

    }
    if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {
        gp <- ggplot(plot_data, aes_string(opt$rdim1, opt$rdim2, color=color_var))
    } else {
        gp <- ggplot(plot_data, aes_string(opt$rdim1, opt$rdim2,
                                           color=color_var, shape=opt$shapefactor))
    }

    if(numeric)
    {
        midpoint <- (max(plot_data[[color_var]]) + min(plot_data[[color_var]]))/2
        gp <- gp + scale_color_gradient2(low="black",mid="yellow",high="red",midpoint=midpoint)
    }

    if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {
        gp <- gp + geom_point(size=opt$pointsize, pch=opt$pointpch,alpha=opt$pointalpha)
    } else {
        gp <- gp + geom_point(size=opt$pointsize, alpha=opt$pointalpha)
    }


    if(color_var == cluster_col)
    {
        gp <- gp + geom_text(data=centers, aes(x, y, label=cluster), color="black")
    }


    # increase the size of the points in the legend.

    if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {

        gp <- gp + guides(color = guide_legend(override.aes = list(size=5, pch=16, alpha=1)))
        gp <- gp + guides(shape = guide_legend(override.aes = list(size=5, pch=16, alpha=1)))


    } else {

        gp <- gp + guides(color = guide_legend(override.aes = list(size=5, alpha=1)))
        gp <- gp + guides(shape = guide_legend(override.aes = list(size=5, alpha=1)))


    }




    ## write out a separate legend if we have more than 20 levels.



    plotfilename = paste(opt$method, color_var, sep=".")
    texCaption <- paste(opt$method,"plot colored by",color_var)

    nlevels <- length(unique(plot_data[[color_var]]))

    if(nlevels > 20 && numeric==FALSE)
    {
        legend <- g_legend(gp)
        gp <- gp + theme(legend.position="none")

        save_ggplots(file.path(opt$outdir, plotfilename),
                     gp,
                     width=7,
                     height=7,
                     to_pdf=opt$pdf)


        legendfilename = paste(opt$method, color_var, "legend", sep=".")
        save_ggplots(file.path(opt$outdir, legendfilename),
                     legend,
                     width=7,
                     height=7,
                     to_pdf=opt$pdf)

        legendCaption <- paste(texCaption, "plot legend")

        tex <- paste(tex,
                     getSubsectionTex(texCaption),
                     getFigureTex(plotfilename,texCaption,
                                  plot_dir_var=opt$plotdirvar),
                     getFigureTex(legendfilename,legendCaption,
                                  plot_dir_var=opt$plotdirvar),
                     sep="\n")

    } else {
        save_ggplots(file.path(opt$outdir, plotfilename),
                     gp,
                     width=7,
                     height=5,
                     to_pdf=opt$pdf)

        tex <- paste(tex,
                     getSubsectionTex(texCaption),
                     getFigureTex(plotfilename,texCaption,
                                  plot_dir_var=opt$plotdirvar),
                     sep="\n")
    }



}


print("Writing out latex snippet")
## write out latex snippet

tex_file <- file.path(opt$outdir,
                      paste("plot.rdims",
                            "factor.tex",
                            sep="."))

writeTex(tex_file, tex)
