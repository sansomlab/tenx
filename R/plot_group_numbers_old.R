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
##           --table=grouping.factors.tsv
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
  require(tenxutils),
  require(colormap)
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
      c("--seuratassay"),
      default="RNA",
      help="The seurat assay from which the counts will be retrieved",
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
opt$subgroupfactor <- strsplit(opt$subgroupfactor, ",")[[1]]
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
    message("doing cell counts and fractions plots for group ",group_var)

    if(opt$subgroupfactor=="none" | !(opt$subgroupfactor %in% colnames(data)))
    {
      message("no subgroup found")

        plot_data <- data %>% group_by_(group_var) %>% summarise(count=n()) %>% mutate(proportion= count/sum(count))

        cm_palette <- colormap(colormap = colormaps$portland,
                             nshade = 2, alpha=0.6)

        gp1 <- ggplot(plot_data, aes_string(group_var, "count"))
        gp1 <- gp1 + geom_bar(stat="identity")

        gp2 <- ggplot(plot_data, aes_string(group_var,"proportion"))
        gp2 <- gp2 + geom_bar(stat="identity")


        gp1 <- gp1 + scale_fill_manual(values=cm_palette) #c("seagreen4","bisque2",")
        gp1 <- gp1 + xlab(group_var)
        gp1 <- gp1 + ylab("number of cells")


        gp2 <- gp2 + scale_fill_manual(values=cm_palette) #c("seagreen4","bisque2",")
        gp2 <- gp2 + xlab(group_var)
        gp2 <- gp2 + scale_y_continuous(labels=scales::percent)
        gp2 <- gp2 + ylab("percentage of cells")

      gp1 <- gp1 + theme_light()
      gp2 <- gp2 + theme_light()

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
                    file.path(opt$outdir,paste0(plotfilename,".plotdata.tsv")),
                    col.names=TRUE,
                    row.names=FALSE,
                    quote=FALSE,
                    sep="\t")

        tex <- c(tex,
                 getSubsectionTex(paste("Cell numbers by",group_var)))

        tex <- c(tex,
                 getFigureTex(plotfilename,
                              plot_title,
                              plot_dir_var=opt$plotdirvar),
                 sep="\n")



    } else {
      for ( subgroup in opt$subgroupfactor){
        message("doing cell counts and fractions plots for group ",group_var, " and subgroup ", subgroup)

        plot_data <- data %>%
              group_by_at(.vars=c(subgroup, group_var)) %>%
              summarise(count=n()) %>%
              mutate(proportion= count/sum(count))

        # plot_data <- data %>% group_by_at(.vars=c(subgroup, group_var)) %>% summarise(count=n()) %>%
        #   group_by_at(.vars=group_var) %>%
        #   mutate(proportion= count/sum(count)) (this would do the percentage for each class within its x variable group)

        # filler<-gtools::mixedsort(unique(as.character(pull(plot_data,subgroup))))
        nsubgroup <- length(unique(plot_data[[subgroup]]))
        cm_palette <- colormap(colormap = colormaps$portland,
                               nshade = nsubgroup, alpha=0.6)


        gp1 <- ggplot(plot_data, aes_string(group_var, "count", group="count", fill=subgroup))
        gp1 <- gp1 + geom_bar(stat="identity", position="dodge")

        gp2 <- ggplot(plot_data, aes_string(group_var, "proportion", group="count", fill=subgroup))
        gp2 <- gp2 + geom_bar(stat="identity", position="dodge")


        gp1 <- gp1 + scale_fill_manual(values=cm_palette)

        gp1 <- gp1 + xlab(group_var)
        gp1 <- gp1 + ylab("number of cells")


        gp2 <- gp2 + scale_fill_manual(values=cm_palette) #scale_fill_viridis_d(breaks=filler, limits=filler)
        gp2 <- gp2 + xlab(group_var)
        gp2 <- gp2 + scale_y_continuous(labels=scales::percent)
        gp2 <- gp2 + ylab("percentage of cells")

        if(nsubgroup > 6 ) {
          gp2 <-gp2 + theme(legend.position = "bottom",
                            legend.key.height = unit(0.1,"in"),
                            legend.text = element_text(size=7))
          gp1 <-gp1 + theme(legend.position = "bottom",
                            legend.key.height = unit(0.1,"in"),
                            legend.text = element_text(size=7))
        }

        gp1 <- gp1 + theme_light()
        gp2 <- gp2 + theme_light()

        gps <- list(gp1, gp2)

        g <- grid.arrange(grobs=gps, ncols=1)

        plot_title <- paste("numbers","group_stats",group_var,"split_by",subgroup,sep=".")

        plotfilename <- paste(plot_title, sep=".") #not necessary

        save_ggplots(file.path(opt$outdir, plotfilename),
                     g,
                     width=7,
                     height=5)

        # save the plot data.
        write.table(plot_data,
                    file.path(opt$outdir,paste0(plotfilename,".plotdata.tsv")),
                    col.names=TRUE,
                    row.names=FALSE,
                    quote=FALSE,
                    sep="\t")

        tex <- c(tex,
                 getSubsectionTex(paste("Cell numbers by",group_var)))

        tex <- c(tex,
                 getFigureTex(plotfilename,
                              plot_title,
                              plot_dir_var=opt$plotdirvar),
                 sep="\n")



        }

    }


    ## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
    ## The palette with grey:

}


# Plot numbers of genes/cell and counts/cell ----

s <- readRDS(opt$seuratobject)

## set the default assay
message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("plot_group_numbers.R running with default assay: ", DefaultAssay(s))

cdata <- GetAssayData(object = s, slot = "counts")

# ngenes <- apply(cdata,2,function(x) sum(x>0))
# ncounts <- apply(cdata,2,sum)
ngenes <- Matrix::colSums(cdata>0)
ncounts <- Matrix::colSums(cdata)


cinfo <- data.frame(ngenes=ngenes,ncounts=ncounts)

data <- merge(data,cinfo,by=0)

## Make one whole-page plot per grouping variable
for(group_var in group_vars)
{
    message("doing genes and counts per-cell plots for group ",group_var)

    if(opt$subgroupfactor=="none" | !(opt$subgroupfactor %in% colnames(data)))
    {
      message("no subgroup found")

      cm_palette <- colormap(colormap = colormaps$portland,
                               nshade = 2, alpha=0.6)

      gp1 <- ggplot(data, aes_string(group_var,"ngenes"))  + geom_boxplot()
      gp1 <-gp1 + xlab(group_var) + ylab("number of genes per cell")

      gp2 <- ggplot(data, aes_string(group_var,"ncounts"))  + geom_boxplot()
      gp2 <- gp2 + xlab(group_var) + ylab("raw counts per cell")

      gp1 <- gp1 + scale_fill_manual(values=cm_palette) ## c("seagreen4","bisque2"))
      gp2 <- gp2 + scale_fill_manual(values=cm_palette) ## c("seagreen4","bisque2"))

      gp1 <- gp1 + theme_light()
      gp2 <- gp2 + theme_light()


      gps <- list(gp1, gp2)

      g <- grid.arrange(grobs=gps, ncols=1)

      plot_title <- paste("numbers", "cell_stats", group_var, sep=".")

      plotfilename <- plot_title

      save_ggplots(file.path(opt$outdir, plotfilename),
                   g,
                   width=6,
                   height=8)

      tex <- c(tex,
               getSubsectionTex(paste("Gene and UMIs by",group_var)))

      tex <- c(tex,
               getFigureTex(plotfilename, plot_title,
                            plot_dir_var=opt$plotdirvar),
               sep="\n")


    } else {
      for ( subgroup in opt$subgroupfactor){
        message("doing cell counts and fractions plots for group ",group_var, " and subgroup ", subgroup)

                                        # filler <- gtools::mixedsort(unique(as.character(data[,subgroup])))

        nsubgroup <- length(unique(plot_data[[subgroup]]))
        cm_palette <- colormap(colormap = colormaps$portland,
                               nshade = nsubgroup, alpha=0.6)

        gp1 <- ggplot(data, aes_string(group_var,"ngenes" ,fill=subgroup))  + geom_boxplot()
        gp1 <-gp1 + xlab(group_var) + ylab("number of genes per cell")

        gp2 <- ggplot(data, aes_string(group_var,"ncounts" ,fill=subgroup))  + geom_boxplot()
        gp2 <- gp2 + xlab(group_var) + ylab("raw counts per cell")
        gp1 <- gp1 + scale_fill_manual(values=cm_palette)
        gp2 <- gp2 + scale_fill_manual(values=cm_palette)
        gp1 <- gp1 + theme_light()
        gp2 <- gp2 + theme_light()

        if(nsubgroup > 6 ) {
          gp2 <-gp2 + theme(legend.position = "bottom",
                            legend.key.height = unit(0.1,"in"),
                            legend.text = element_text(size=7))
          gp1 <-gp1 + theme(legend.position = "bottom",
                            legend.key.height = unit(0.1,"in"),
                            legend.text = element_text(size=7))
        }

        gp1 <- gp1 + theme_light()
        gp2 <- gp2 + theme_light()

        gps <- list(gp1, gp2)

        g <- grid.arrange(grobs=gps, ncols=1)

        plot_title <- paste("numbers", "cell_stats", group_var, "split_by", subgroup, sep=".")

        plotfilename <- plot_title

        save_ggplots(file.path(opt$outdir, plotfilename),
                     g,
                     width=6,
                     height=8)

        tex <- c(tex,
                 getSubsectionTex(paste("Gene and UMIs by",group_var)))

        tex <- c(tex,
                 getFigureTex(plotfilename, plot_title,
                              plot_dir_var=opt$plotdirvar),
                 sep="\n")

      }

    }


}

## write out latex snippet
writeTex(file.path(opt$outdir, paste("number","plots","tex", sep=".")),
         tex)

message("plot_group_numbers.R final default assay: ", DefaultAssay(s))

message("completed")
