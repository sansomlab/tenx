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

    if(length(intersect(rownames(plot_data), rownames(meta_data))) < length(rownames(meta_data)))
    {
        stop("Not all cell barcodes are present in the given metadata table")
    } else {

        meta_data <- meta_data[rownames(plot_data),]
        plot_data <- cbind(plot_data, meta_data)
    }
}

if ("cluster" %in% colnames(plot_data)){
  plot_data$cluster <- factor(plot_data$cluster, levels=sort(as.numeric(unique(plot_data$cluster))))
}


color_vars <- strsplit(opt$colorfactors,",")[[1]]
tex = ""

print("Making tSNE plots colored by each of the factor variables")
## Make one whole-page tSNE plot per color variable
for(color_var in color_vars)
{
    print(paste("Making",color_var,"tSNE plot"))

    ## If a variable comprises only integers, coerce it to a character vector

    numeric = FALSE
    if(is.numeric(plot_data[[color_var]]))
    {
        if(all(plot_data[[color_var]] == round(plot_data[[color_var]]))
           & length(unique(plot_data[[color_var]])) < 50)
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

    gp <- gp + geom_point(size=opt$pointsize)

    plotfilename = paste(opt$method, color_var, sep=".")

    save_ggplots(file.path(opt$outdir, plotfilename),
                 gp,
                 width=6,
                 height=4,
                 to_pdf=opt$pdf)

    texCaption <- paste(opt$method,"plot colored by",color_var)

    tex <- paste(tex,
                 getSubsectionTex(texCaption),
                 getFigureTex(plotfilename,texCaption,
                              plot_dir_var=opt$plotdirvar),
                 sep="\n")

}


print("Writing out latex snippet")
## write out latex snippet

tex_file <- file.path(opt$outdir,
                      paste("plot.rdims",
                            "factor.tex",
                            sep="."))

writeTex(tex_file, tex)
