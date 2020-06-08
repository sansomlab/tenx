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
  require(Seurat),
  require(tenxutils),
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
      c("--seuratobject"),
      default="none",
      help="The seurat object (e.g. begin.rds)"
    ),
    make_option(
      c("--seuratassay"),
      default="RNA",
      help="The seurat assay to pull the expression data from"
    ),
    make_option(
      c("--plotdirvar"),
      default="clusterMarkerTSNEPlotsDir",
      help="The name of the latex var specifying the location of the plots"
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
      c("--outdir"),
      default="seurat.out.dir",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


## read in the raw count matrix
s <- readRDS(opt$seuratobject)

## set the default assay
message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("plot_rdims_gene.R running with default assay: ", DefaultAssay(s))


data <- GetAssayData(object = s, slot = "data")

if("gene_id" %in% colnames(s@misc))
{
    ## relabel data with ensembl ids
    rownames(data) <- s@misc[rownames(data),"gene_id"]
    id_col = "gene_id"
} else {
    ## if "gene_id" not avaliable fall back to "gene"
    id_col = "gene"
}

##
plot_data <- read.table(opt$table,sep="\t",header=TRUE)

print(head(plot_data))
rownames(plot_data) <- plot_data$barcode

## read in the table containing the genes to visualise
genes <- read.table(opt$genetable, header=TRUE, sep="\t",as.is=TRUE)
if("common_name" %in% colnames(genes))
{
    plot_names <- make.unique(paste(genes$common_name," (",genes$gene_name,")",sep=""))
} else {
    plot_names <- make.unique(genes$gene_name)
    }
genes$plot_name <- plot_names

rownames(genes) <- genes[[id_col]]

## get the geneset name
geneset <- gsub(".txt", "", basename(opt$genetable))
geneset <- gsub(".csv", "", geneset)

## initialise the tex snippet
tex = "" # getSubsectionTex(geneset)

good_ids <- genes[[id_col]][genes[[id_col]] %in% rownames(data)]

## subset the gene expression matrix
exprs <- as.matrix(data[good_ids, plot_data$barcode, drop=F])

print(dim(exprs))

rownames(exprs) <-genes$plot_name[match(rownames(exprs),genes[[id_col]])]

## get the 95th quantile of log10 n+1 counts > 0
## upper_limit = quantile(log10(data[data>0]+1),0.95)

# upper_limit = quantile(data[data>0],0.95)
# logical crashes transforming sparse matrix into dense. Use indexes, load Matrix.

xx<-which(data>0)
upper_limit = quantile(data[xx], 0.95)

## log transform
## exprs <- log10(exprs+1)
exprs[exprs > upper_limit] <- upper_limit

plot_genes <- rownames(exprs)

ngenes = length(plot_genes)
npages = ceiling(ngenes/9)

message("processing ", npages," pages (", ngenes, " genes)")


seqWrapper <- function(lb, ub, by=1) {
    s <- c()
    if(!ub < lb) s <- seq(lb,ub, by=by)
    return(s)
}

for(page in seqWrapper(1,npages))
{
    message("working on page: ", page)
    start <- ((page-1)*9) + 1
    end <- min((page*9),ngenes)
    nplots = (end - start) + 1

    message("number of plots to make: ", nplots)
    nplotrows = ceiling(nplots/3)
    genes_of_interest <- plot_genes[start:end]

    # subset the table with the reduced coordinates
    if(tolower(opt$shapefactor)!="none"){
        plot_data <- plot_data[,c(opt$rdim1,opt$rdim2,opt$shapefactor)]
    } else {
        plot_data <- plot_data[, c(opt$rdim1,opt$rdim2)]
    }

    # add the gene expression information
    plot_data <- merge(plot_data, t(exprs[genes_of_interest,]), by=0)
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL

    # melt the data for plotting by gene
    if(tolower(opt$shapefactor)!="none"){
        melted_plot_data <- melt(plot_data, id.vars=c(opt$rdim1,opt$rdim2,opt$shapefactor))
    } else {
        melted_plot_data <- melt(plot_data, id.vars=c(opt$rdim1,opt$rdim2))
    }

    theme_set(theme_gray(base_size = 8))
    if(tolower(opt$shapefactor)=="none" | !opt$shapefactor %in% colnames(melted_plot_data))
    {
        gp <- ggplot(melted_plot_data, aes_string(opt$rdim1, opt$rdim2,
                                                  color="value"))
    } else {
        gp <- ggplot(melted_plot_data, aes_string(opt$rdim1, opt$rdim2,
                                                  color="value",
                                                  shape=opt$shapefactor))
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

    gp <- gp + geom_point(size=opt$pointsize)
    gp <- gp + theme(panel.background = element_rect(fill = 'lightgrey'))

    if(nplotrows==1) { h=3.5 }
    if(nplotrows==2) { h=8 }
    if(nplotrows==3) { h=10 }

    pagename <- paste("genes",page,sep="_")
    plotfilename <- paste("plot.rdims.genes", geneset, page, sep=".")

    save_ggplots(file.path(opt$outdir, plotfilename),
                 gp,
                 width=10,
                 height=h,
                 dpi=300)

    texCaption <- paste(geneset," genes (",page," of ",npages,")",sep="")
    tex <- c(tex, getFigureTex(plotfilename, texCaption,
                               plot_dir_var=opt$plotdirvar))
}

message("Completed plotting.")

tex_file <- file.path(opt$outdir, paste("plot.rdims.genes", geneset, "tex", sep="."))

writeTex(tex_file, tex)

message("plot_rdims_gene.R final default assay: ", DefaultAssay(s))

message("completed")
