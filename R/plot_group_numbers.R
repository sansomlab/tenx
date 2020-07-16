## Title ----
##
## Visualise the numbers of single-cells falling into different groups
##
## Description ----
##

# Libraries ----

stopifnot(
  require(optparse),
  require(ggplot2),
  require(reshape2),
  require(dplyr),
  require(tenxutils),
  require(colormap)
)

# Options ----

option_list <- list(
    make_option(
      c("--metadata"),
      default="none",
      help="A table containing the metadata with grouping information"),
    make_option(
      c("--clusters"),
      default="none",
      help="A table containing the clusters"),
    make_option(
      c("--seuratobject"),
      default="none",
      help="The seurat object (typically begin.rds)"
    ),
    make_option(
      c("--seuratassay"),
      default="RNA",
      help="The seurat assay from which the counts will be retrieved",
      ),
    make_option(
      c("--title"),
      default="none",
      help="Used for the filename"
      ),
    make_option(
      c("--group"),
      default="none",
      help="A column in the cell metadata used to define the groups"
      ),
    make_option(
      c("--subgroup"),
      default="none",
      help="Column(s) in the cell metadata for subgroup"
    ),
    make_option(
      c("--subgrouplevels"),
      default=NULL,
      help="the subgroup levels in order, comma separated"
    ),
    make_option(
        c("--facet"),
        default=FALSE,
        action="store_true",
        help="Facet the plot by group."
    ),
    make_option(
        c("--freey"),
        default=FALSE,
        action="store_true",
        help='should the y scale be free when facetting?'
    ),
    make_option(
        c("--ncol"),
        default=5,
        type="integer",
        help='should the y scale be free when facetting?'
    ),
    make_option(
        c("--replicate"),
        default="none",
        help="Column in the cell metadata used to define the replicates"
    ),
    make_option(
        c("--xlab"),
        default=NULL,
        help="x axis label"
    ),
    make_option(
        c("--ylab"),
        default=NULL,
        help="y axis label"
    ),
    make_option(
        c("--stat"),
        default="none",
        help=paste('The statistic to be plotted. Either:',
                   '(1) A column name in the metadata',
                   '(2) "total_UMI" or "ngenes" which will be computed from the seurat object, if not in the metadata',
                   '(3) the summary stats "count" or "pct" which will be computed by dplyr'),
    ),
    make_option(
      c("--geom"),
      default="bar",
      help='"bar" or "boxplot"'
    ),
    make_option(
      c("--trans"),
      default=NULL,
      help='Currently only sqrt is supported, applied to y'
      ),
    make_option(
      c("--width"),
      default=7,
      type="double",
      help="the width of the plot"
      ),
    make_option(
      c("--height"),
      default=7,
      type="double",
      help="the width of the plot"
      ),
    make_option(
      c("--plotdirvar"),
      default="groupNumbersDir",
      help="the latex var holding the name of the dir with the plots"
      ),
    make_option(
      c("--outdir"),
      default="seurat.out.dir",
      help="outdir"),
    make_option(
      c("--pdf"),
      default=FALSE,
      help="outdir"),
    make_option(
      c("--savedata"),
      default=TRUE,
      help="should the plot data be saved?")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

## Read in the table with the grouping information
data <- read.table(opt$metadata,sep="\t",header=T)
rownames(data) <- data$barcode

if(opt$cluster!="none")
{


    clusters <- read.table(opt$clusters, sep="\t", header=TRUE)
    rownames(clusters) <- clusters$barcode
    clusters$barcode <- NULL

    if(!all(rownames(data) %in% rownames(clusters)))
    {
        stop("Not all cell barcodes are present in the given cluster table")
    }

    data$cluster <- clusters[rownames(data), "cluster_id"]


}




## ########################################################################### #
## ########################## Add missing stats ############################## #
## ########################################################################### #

if(opt$stat %in% c("total_UMI", "ngenes"))
{
    if(!opt$stat %in% colnames(data))
    {
        require(Seurat)
        # add the statistic from the seurat object
        if (endsWith(opt$seuratobject, ".rds")) {
          message(sprintf("readRDS: %s", opt$seuratobject))
          s <- readRDS(opt$seuratobject)
        } else {
          require(SeuratDisk)
          message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
          s <- LoadH5Seurat(opt$seuratobject)
        }
      
        ## set the default assay
        message("Setting default assay to: ", opt$seuratassay)
        DefaultAssay(s) <- opt$seuratassay

        message("plot_group_numbers.R running with default assay: ", DefaultAssay(s))

        cdata <- GetAssayData(object = s, slot = "counts")

        if( opt$stat=="ngenes")
        {
            ngenes <- Matrix::colSums(cdata>0)
            data$ngenes <- ngenes[rownames(data)]
        }

        else if(opt$stat=="total_UMI")
        {
            total_UMI <- Matrix::colSums(cdata)
            data$total_UMI <- total_UMI[rownames(data)]
        } else { stop("stat not recognised") }
        # free the memory.
        rm(s)
    }
}

## ########################################################################### #
## ################ summarise the data for plotting ########################## #
## ########################################################################### #

# define the grouping variables
group_vars = c(opt$group, opt$subgroup, opt$replicate)
group_vars <- group_vars[!group_vars=="none"]

message("using group variables")
print(group_vars)

## check the group variables are in the data columns
for(gvar in group_vars)
{
    if(!gvar %in% colnames(data))
    {
        stop(paste("Grouping variable: ", gvar," not found in data"))
    }
}


if(opt$stat=="count")
{
    plot_data <- data %>%
        group_by_at(group_vars) %>%
        summarise(count=n())

} else if (opt$stat == "pct") {
    if(opt$replicate=="none")
    {
        ## this the within group percent.
        ## to normalise across groups, specify a replicate factor.

        plot_data <- data %>%
            group_by_at(group_vars) %>%
            summarise(count=n()) %>%
            mutate(pct=(count/sum(count))*100)

    } else {
        plot_data <- data %>% add_count((!!as.name(opt$replicate)), name="rep_n")
        message("head plot data after adding replicate count column")

        print(plot_data[,c(colnames(plot_data)[1:5],"rep_n")])

        plot_data <- plot_data %>%
            group_by_at(group_vars) %>%
            summarise(count=n(), n=mean(rep_n)) %>%
            mutate(pct=(count/n)*100)

        ## check that the percentages by replicate add to 100
        check <- plot_data %>%
            group_by_at(opt$replicate) %>%
            summarise(total_pct=sum(pct))

        if(mean(check$total_pct) != 100) { stop("pcts do not sum to 100") }

        message("head of pct check tibble")
        print(check)

    }
} else if (!opt$stat %in% colnames(data)) {
    stop("The requested statistic is not recognised or present in the metadata")
} else {

    if(opt$replicate=="none")
    {
        stop("A replicate column must be specified when plotting by a named column")

        } else {
    ## no need to summarise, plotting boxplots on column in data.

            plot_data <- data
            }
}

message("head of the plot data..")
print(head(plot_data))


## ########################################################################### #
## ######################## draw out the plots ############################### #
## ########################################################################### #


## set the x and fill variables.
if(opt$subgroup == "none")
{
            subgroup <- FALSE
            xvar <- opt$group
            fvar <- opt$group

} else {
    subgroup <- TRUE

    if(!is.null(opt$subgrouplevels))
        {
            sgl = strsplit(gsub(" ","", opt$subgrouplevels),",")[[1]]

            if(!all(sgl  %in% unique(plot_data[[opt$subgroup]])))
            {
                stop("given levels do not match levels present in subgroup")
            }

            plot_data[[opt$subgroup]] <- factor(plot_data[[opt$subgroup]],
                                                levels = sgl)
        }

    fvar <- opt$subgroup
    if(opt$facet)
        {
            xvar <- opt$subgroup
        } else { xvar <- opt$group }

}

message("fvar: ", fvar)
message("xvar: ", xvar)

## convert numeric x factor to char with correct ordering
if(is.numeric(plot_data[[xvar]]))
{

    x <- plot_data[[xvar]]
    num_levels = unique(x)
    char_levels = as.character(num_levels)
    char_levels = char_levels[order(num_levels)]

    plot_data[[xvar]] <- factor(as.character(plot_data[[xvar]]),
                                   levels=char_levels)

}

gp <- ggplot(plot_data, aes_string(xvar, opt$stat, fill=fvar))

if(opt$geom == "bar")
{
    message(" drawing bar plots")

    gp <- gp + geom_bar(stat="identity", position="dodge")


} else if(opt$geom == "boxplot") {
    message("drawing box plots")

    gp <- gp + geom_boxplot()

} else { stop("geom not supported") }

## optional facetting
if(subgroup && opt$facet) {
    if(opt$freey)
        {
            gp <- gp + facet_wrap(formula(paste("~",opt$group,sep="")), scales="free_y", ncol=opt$ncol)
        } else {
            gp <- gp + facet_wrap(formula(paste("~",opt$group,sep="")),ncol=opt$ncol)
            }

}

gp <- gp + theme_light()

## sort out the colors
nsubgroup <- length(unique(plot_data[[fvar]]))
cm_palette <- colormap(colormap = colormaps$portland,
                           nshade = nsubgroup, alpha=0.6)

gp <- gp + scale_fill_manual(values=cm_palette) #c("seagreen4","bisque2",")

## optional y axis transformation
if(!is.null(opt$trans))
{
    if(opt$trans=="sqrt")
    {    gp <- gp + scale_y_continuous(trans=opt$trans)
    } else { stop("given transformation not implemented") }
}


## tidy up the axis labels.
gp <- gp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

if(!is.null(opt$xlab)) { gp <- gp + xlab(opt$xlab) }
if(!is.null(opt$ylab)) { gp <- gp + ylab(opt$ylab) }


## save the plots
save_ggplots(file.path(opt$outdir, opt$title),
             gp,
             width=opt$width,
             height=opt$height,
             to_pdf=opt$pdf)


## optionally save the plot data
if(opt$savedata)
{
    data_file <- file.path(opt$outdir,
                           paste(opt$title, "data.tsv.gz", sep="."))
    write.table(plot_data, gzfile(data_file),
                sep="\t", col.names=T, row.names=F, quote=F)
}

message("completed")
