## Title ----
##
## Visualise the numbers of single-cells falling into different groups
##
## Description ----
##
## For the given grouping variables, make plots of
## (i) numbers of cells
## (ii) fractions of cells
## (iii) ngenes/cell
## (iv) ncounts/cell.
##
## Details ----
##
## Only grouping factors with more than one level are plotted.
##
## Usage ----
##
## $ Rscript getGenesetAnnotations.R
##           --table=grouping.factors.txt
##           --seuratobject=seurat.Robj
##           --groupfactors=cluster
##           --subgroupfactors=patient
##           --outdir=.

# Libraries ----

stopifnot(
  require(optparse),
  require(ggplot2),
  require(reshape2),
  require(dplyr),
  require(Seurat),
  require(gridExtra),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
      c("--table"),
      default="none",
      help="A table containing the grouping information"),
    make_option(
      c("--seuratobject"),
      default="none",
      help="The seurat object (typically begin.rds)"
      ),
    make_option(
      c("--subgroupfactor"),
      default="none",
      help="A column in the cell metadata used to define subgroups e.g. the condition"
      ),
    make_option(
      c("--groupfactors"),
      default="none",
      help="Column(s) in the cell metadata to use for grouping, e.g. cluster or State"
    ),
    make_option(
      c("--plotdirvar"),
      default="groupNumbersDir",
      help="the latex var holding the name of the dir with the plots"
      ),
    make_option(
      c("--outdir"),
      default="seurat.out.dir",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

## Read in the table with the grouping information
data <- read.table(opt$table,sep="\t",header=T)
rownames(data) <- data$barcode

group_vars_all <- strsplit(opt$groupfactors,",")[[1]]
tex = ""

## Drop grouping variables
## that only have one level
group_vars <- c()
for(group_var in group_vars_all)
{
    data[[group_var]] <- as.factor(data[[group_var]])

    if(length(levels(data[[group_var]]))>1)
    {
        group_vars <- c(group_vars, group_var)
    }
}

# Plot cell counts and fractions ----
## Make one whole-page plot per grouping variable
for(group_var in group_vars)
{

    if(opt$subgroupfactor=="none" | !(opt$subgroupfactor %in% colnames(data)))
    {

        plot_data <- data %>% group_by_(group_var) %>% summarise(count=n()) %>% mutate(proportion= count/sum(count))

        gp1 <- ggplot(plot_data, aes_string(group_var, "count"))
        gp1 <- gp1 + geom_bar(stat="identity")

        gp2 <- ggplot(plot_data, aes_string(group_var,"proportion"))
        gp2 <- gp2 + geom_bar(stat="identity")


    } else {

        plot_data <- data %>% group_by_(opt$subgroupfactor, group_var) %>% summarise(count=n()) %>% mutate(proportion= count/sum(count))


        gp1 <- ggplot(plot_data, aes_string(group_var, "count", group="count", fill=opt$subgroupfactor))
        gp1 <- gp1 + geom_bar(stat="identity", position="dodge")

        gp2 <- ggplot(plot_data, aes_string(group_var, "proportion", group="count", fill=opt$subgroupfactor))
        gp2 <- gp2 + geom_bar(stat="identity", position="dodge")

    }


    ## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    ## The palette with grey:

    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    gp1 <- gp1 + scale_fill_manual(values=cbPalette) #c("seagreen4","bisque2",")
    gp1 <- gp1 + xlab(group_var)
    gp1 <- gp1 + ylab("number of cells")


    gp2 <- gp2 + scale_fill_manual(values=cbPalette) #c("seagreen4","bisque2",")
    gp2 <- gp2 + xlab(group_var)
    gp2 <- gp2 + scale_y_continuous(labels=scales::percent)
    gp2 <- gp2 + ylab("percentage of cells")

    gps <- list(gp1, gp2)

    g <- grid.arrange(grobs=gps, ncols=1)

    plot_title <- paste("numbers","group_stats",group_var,sep=".")

    plotfilename <- paste(plot_title, sep=".")

    save_ggplots(file.path(opt$outdir, plotfilename),
                 g,
                 width=6,
                 height=4)

    # save the plot data.
    write.table(plot_data,
                file.path(opt$outdir,paste0(plotfilename,".plotdata.txt")),
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE,
                sep="\t")

    tex <- paste(tex,
                 getFigureTex(plotfilename,
                              plot_title,
                              plot_dir_var=opt$plotdirvar),
                 sep="\n")

}


# Plot numbers of genes/cell and counts/cell ----

s <- readRDS(opt$seuratobject)
cdata <- s@raw.data

ngenes <- apply(cdata,2,function(x) sum(x>0))
ncounts <- apply(cdata,2,sum)
cinfo <- data.frame(ngenes=ngenes,ncounts=ncounts)

data <- merge(data,cinfo,by=0)

## Make one whole-page plot per grouping variable
for(group_var in group_vars)
{

    if(opt$subgroupfactor=="none" | !(opt$subgroupfactor %in% colnames(data)))
    {

        gp1 <- ggplot(data, aes_string(group_var,"ngenes"))  + geom_boxplot()
        gp1 <-gp1 + xlab(group_var) + ylab("number of genes per cell")

        gp2 <- ggplot(data, aes_string(group_var,"ncounts"))  + geom_boxplot()
        gp2 <- gp2 + xlab(group_var) + ylab("raw counts per cell")

    } else {
        gp1 <- ggplot(data, aes_string(group_var,"ngenes" ,fill=opt$subgroupfactor))  + geom_boxplot()
        gp1 <-gp1 + xlab(group_var) + ylab("number of genes per cell")

        gp2 <- ggplot(data, aes_string(group_var,"ncounts" ,fill=opt$subgroupfactor))  + geom_boxplot()
        gp2 <- gp2 + xlab(group_var) + ylab("raw counts per cell")

    }

    gp1 <- gp1 + scale_fill_manual(values=cbPalette) ## c("seagreen4","bisque2"))
    gp2 <- gp2 + scale_fill_manual(values=cbPalette) ## c("seagreen4","bisque2"))

    gps <- list(gp1, gp2)

    g <- grid.arrange(grobs=gps, ncols=1)

    plot_title <- paste("numbers", "cell_stats", group_var, sep=".")

    plotfilename <- plot_title

    save_ggplots(file.path(opt$outdir, plotfilename),
                 g,
                 width=6,
                 height=8)

    tex <- paste(tex,
                 getFigureTex(plotfilename, plot_title,
                              plot_dir_var=opt$plotdirvar),
                 sep="\n")

}

## write out latex snippet
writeTex(file.path(opt$outdir, paste("number","plots","tex", sep=".")),
         tex)
