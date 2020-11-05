## Title ----
##
## Visualise gene expression levels on reduced dimensions plots
##
## Description ----
##
## This script visualises per-cell gene expression levels on a 2D
## reduced dimensions plot.
##
## Details ----
##
## Expression levels are taken from the "data" slot,
## and truncated at the 95th expression percentile
##
## Usage ----
##
## See options.

# Libraries ----

stopifnot(
  require(optparse),
  require(ggplot2),
  require(reshape2),
  require(tenxutils),
  require(BiocParallel),
  require(Matrix)
)

# Options ----

option_list <- list(
    make_option(
      c("--table"),
      default="none",
      help="A table containing the reduced coordinates and phenotype information"
      ),
    make_option(
      c("--method"),
      default="tSNE",
      help="Normally the type of dimension reduction"
      ),
    make_option(
      c("--rdim1"), default="tSNE_1",
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
      c("--genetable"),
      default="none",
      help="A tab-delimited text file containing gene_id (or gene if s@misc$gene_id is not set) and gene_name"
    ),
    make_option(
      c("--cluster"),
      default="none",
      help="the cluster for which to make the plots"
      ),
    make_option(
      c("--name"),
      default="none",
      help="the name of the set of genes"
      ),
    make_option(
      c("--assaydata"),
      default="none",
      help="The seurat object (e.g. begin.rds)"
    ),
    make_option(
      c("--pointsize"),
      default=0.5,
      help="The point size for the tSNE plots"
    ),
    make_option(
      c("--pointpch"),
      default="16",
      help="The point shape if no shape factor is set for the tSNE plots"
      ),

    make_option(
      c("--pointalpha"),
      default=0.8,
      help="The alpha setting for the points on the tSNE plots"
    ),
    make_option(
      c("--pdf"),
      default=FALSE,
      help="produce pdf plots"
    ),
    make_option(
        c("--maxcells"),
        default=Inf,
        type="integer",
        help="max cells to plot, if there are more they will be randomly downsampled"),
    make_option(
      c("--outdir"),
      default="seurat.out.dir",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

if(opt$pointpch != ".") { opt$pointpch <- as.numeric(opt$pointpch) }

cat("Running with options:\n")
print(opt)


## read in the assay data
#data <- read.table(opt$assaydata, sep="\t", header=T, as.is=T, check.names=F)
message("reading in assay data")
data <- loadSeurat(opt$assaydata)

message("reading in the coordinates")
plot_data <- read.table(opt$table,sep="\t",header=TRUE)

rownames(plot_data) <- plot_data$barcode

if(nrow(plot_data) > opt$maxcells)
{
    ## randomly downsample the plot data.
    message("downsampling the plot data")
    plot_data <- plot_data[sample(rownames(plot_data), opt$maxcells, replace=FALSE),]
}

message("reading in the features")
## read in the table containing the genes to visualise
genes <- read.table(opt$genetable, header=TRUE, sep="\t",as.is=TRUE)
rownames(genes) <- genes$gene

genes <- genes[genes$cluster==opt$cluster,]

good_ids <- genes$gene[genes$gene %in% rownames(data)]


#print(head(data))
data <- as.matrix(data)
xx<-which(data>0)
upper_limit = quantile(data[xx], 0.95)

## subset the gene expression matrix
exprs <- as.matrix(data[good_ids, plot_data$barcode, drop=F])
print(dim(exprs))

rm(data)

exprs[exprs > upper_limit] <- upper_limit

plot_genes <- rownames(exprs)

ngenes = length(plot_genes)
npages = ceiling(ngenes/9)

message("processing ", npages," pages (", ngenes, " genes)")

options(bitmapType='cairo')

seqWrapper <- function(lb, ub, by=1) {
    s <- c()
    if(!ub < lb) s <- seq(lb,ub, by=by)
    return(s)
}


message("subset table")

## subset the table with the reduced coordinates
if(tolower(opt$shapefactor)!="none"){
    plot_data <- plot_data[,c(opt$rdim1,opt$rdim2,opt$shapefactor)]
} else {
    plot_data <- plot_data[, c(opt$rdim1,opt$rdim2)]
}


draw_page <- function(page,
                      shapefactor=opt$shapefactor,
                      ngenes=ngenes,
                      plot_genes=plot_genes,
                      exprs=exprs,
                      pointpch=opt$pointpch,
                      pointsize=opt$pointsize,
                      pointalpha=opt$pointalpha,
                      rdim1 = opt$rdim1,
                      rdim2 = opt$rdim2,
                      outdir = opt$outdir,
                      to_pdf = opt$pdf
                      )
{

    message(paste("working on page: > ", page, " < "))
    start <- ((page-1)*9) + 1
    end <- min((page*9),ngenes)
    nplots = (end - start) + 1

    message(paste("number of plots to make: > ", nplots, " < "))
    nplotrows = ceiling(nplots/3)
    genes_of_interest <- plot_genes[start:end]


    message("add gene information")
    # add the gene expression information
    plot_data <- merge(plot_data, t(exprs[genes_of_interest,]), by=0)
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL


    message("melt the data")
    # melt the data for plotting by gene
    if(tolower(shapefactor)!="none"){
        melted_plot_data <- melt(plot_data, id.vars=c(rdim1,rdim2,shapefactor))
    } else {
        melted_plot_data <- melt(plot_data, id.vars=c(rdim1,rdim2))
    }

    message("draw the plot")
    theme_set(theme_gray(base_size = 8))
    if(tolower(shapefactor)=="none" | !shapefactor %in% colnames(melted_plot_data))
    {
        gp <- ggplot(melted_plot_data, aes_string(rdim1, rdim2,
                                                  color="value"))
    } else {
        gp <- ggplot(melted_plot_data, aes_string(rdim1, rdim2,
                                                  color="value",
                                                  shape=shapefactor))
    }

    if(length(unique(melted_plot_data$variable)) > 1)
        {
            gp <- gp + facet_wrap(~variable, scales="free", ncol=3)
        }

    gp <- gp + scale_color_gradientn(colors=c("#ffffcc", "#ffeda0", "#fed976", "#feb24c",
                                              "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026",
                                              "#800026"),
                                     name="log\nnormalised\ncounts+1)",
                                     limits=c(0,upper_limit))


    if(tolower(shapefactor)=="none" | !shapefactor %in% colnames(melted_plot_data))
    {
        gp <- gp + geom_point(size=pointsize, pch=pointpch, alpha=pointalpha)
    } else {
        gp <- gp + geom_point(size=pointsize, alpha=pointalpha)
    }

    gp <- gp + theme(panel.background = element_rect(fill = 'lightgrey'))

    if(nplotrows==1) { h=3.5 }
    if(nplotrows==2) { h=8 }
    if(nplotrows==3) { h=10 }

    pagename <- paste("genes",page,sep="_")
    plotfilename <- paste(opt$name, "cluster", opt$cluster,
                          "page", page, sep=".")

    message("save the plot")
    save_ggplots(file.path(outdir, plotfilename),
                 gp,
                 width=10,
                 height=h,
                 dpi=300,
                 to_pdf=to_pdf)

    return(page)
}

xpages <- seqWrapper(1,npages)
print(xpages)

message("drawing the pages")

for(page in xpages)
{
    draw_page(page,
              shapefactor=opt$shapefactor,
              ngenes=ngenes,
              plot_genes=plot_genes,
              exprs=exprs,
              pointpch=opt$pointpch,
              pointsize=opt$pointsize,
              pointalpha=opt$pointalpha,
              rdim1=opt$rdim1,
              rdim2=opt$rdim2,
              outdir=opt$outdir,
              to_pdf=opt$pdf)
}

message("Completed plotting.")
message("completed")
