##
## Make a summary plot of the numbers of DE genes per cluster
##

# Libraries ----

stopifnot(
  require(optparse),
  require(gplots),
  require(reshape2),
  require(xtable),
  require(ggplot2),
  require(ggrepel),
  require(gridExtra),
  require(dplyr),
  require(tenxutils)
)

# Options ----

## run only specified comparison (in order for parallel execution)
## deal with the options #mindiffpct, minpct, thresh.use
option_list <- list(
    make_option(c("--degenes"), default="begin.Robj",
                help="Summary table of differentially expressed genes"),
    make_option(c("--clusterids"), default="begin.Robj",
                help="clusterids"),
    make_option(c("--testfactor"), default="NULL",
                help="a metadata factor used to group the violin plots"),
    make_option(c("--a"), default="NULL",
                help="first level of constrast"),
    make_option(c("--b"), default="NULL",
                help="second level of constrast"),
    make_option(c("--minfc"), default=2,
                help="Fold change threshold for DE genes (absolute)"),
    make_option(c("--minpadj"), default=0.05,
                help="The threshold for the adjusted p-value"),
    make_option(c("--plotdirvar"), default="clusterMarkerDEPlotsDir",
                help="latex var containing the locaiton of the plots"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )


opt <- parse_args(OptionParser(option_list=option_list))

if(opt$testfactor=="NULL") { opt$testfactor <- NULL }

cat("Running with options:\n")
print(opt)

# set the run specs
run_specs <- paste(opt$numpcs,opt$resolution,opt$algorithm,opt$testuse,sep="_")


## set up the ggplot theme
theme_set(theme_classic(base_size = 10))

if(is.null(opt$testfactor))
{
    file_suffix = NULL
    a = "cluster_mean"
    b = "other_mean"
} else {
    file_suffix = paste0("between")
    a = paste(opt$a,"mean",sep="_")
    b = paste(opt$b,"mean",sep="_")
}

## initialise the text snippet
tex = ""

message("reading in the DE genes")
## read in the de genes
degenes <- read.table(gzfile(opt$degenes),header=T,as.is=T,sep="\t")

message("summarising and melting the data")

# Summarise and melt the data
summarised_data <- degenes %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      positive=length(which(p.adj < opt$minpadj & avg_logFC >= log(opt$minfc))),
      negative=length(which(p.adj < opt$minpadj & avg_logFC <= -log(opt$minfc))))

melted_data <- melt(summarised_data, id=c("cluster"))

# Set cluster to an ordered factor for plotting
clevels <- unique(melted_data$cluster)
melted_data$cluster <- factor(melted_data$cluster, levels=clevels[order(clevels)])

message("drawing summary barplot")

gp <- ggplot(melted_data, aes(cluster, value, fill=variable))
gp <- gp + geom_bar(stat="identity",position="dodge")
gp <- gp + scale_fill_manual(values=c("bisque2","seagreen4"))
gp <- gp + ylab(paste0("no. genes (p.adj < ", opt$minpadj,
                       ", fold change > ",  opt$minfc, ")"))

nsfn <- paste(c("deNumbers",file_suffix),collapse=".")
nsfp <- file.path(opt$outdir, nsfn)

save_ggplots(nsfp,
             gp,
             width=6,
             height=6)

subsectionTitle <- getSubsubsectionTex(paste("Summary of numbers of DE genes per-cluster"))
tex <- c(tex, subsectionTitle)

deCaption <- paste0("Numbers of differentially expressed genes ",
                    "(adjusted p-value < ", opt$minpadj,
                    ", fold change > ", opt$minfc, ") per cluster")

tex <- c(tex, getFigureTex(nsfn, deCaption,
                           plot_dir_var=opt$plotdirvar))

tex_file <- file.path(opt$outdir,
                      paste(c("characterise.degenes",file_suffix,"tex"),collapse="."))

writeTex(tex_file, tex)
