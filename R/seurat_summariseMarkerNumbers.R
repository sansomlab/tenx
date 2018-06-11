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
    make_option(c("--useminfc"), default=FALSE,
                help="Use minimum fold change for second page of violin plots"),
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

clevels <- unique(degenes$cluster)

## Make a bar plot summarising the number of
## DE genes per cluster
## (rewritten to avoid use of dplyr, which was causing random crashing)

restmp <- c()
for(cluster in clevels)
{
    tmp <- degenes[degenes$p.adj < 0.05 & degenes$cluster==cluster,]

    npos <- nrow(tmp[tmp$avg_logFC > log(1),])
    nneg <- nrow(tmp[tmp$avg_logFC < -log(1),])


    restmp <- c(restmp,
                c(cluster, npos, "positive"),
                c(cluster, nneg, "negative"))
}

nsummary <- data.frame(matrix(restmp,ncol=3,byrow=TRUE))
colnames(nsummary) <- c("cluster","n","type")

message("drawing summary barplot")

gp <- ggplot(nsummary, aes(cluster, n, fill=type))
gp <- gp + geom_bar(stat="identity",position="dodge")
gp <- gp + scale_fill_manual(values=c("seagreen4","bisque2"))
gp <- gp + ylab("no. genes (p.adj < 0.05, fold change > 2)")

nsfn <- paste(c("deNumbers",file_suffix),collapse=".")
nsfp <- file.path(opt$outdir, nsfn)

save_ggplots(nsfp,
             gp,
             width=6,
             height=6)

subsectionTitle <- getSubsubsectionTex(paste("Summary of numbers of DE genes per-cluster"))
tex <- c(tex, subsectionTitle)

deCaption <- "Numbers of differentially expressed genes (adjusted p-value < 0.05, fold change > 2) per cluster"
tex <- c(tex, getFigureTex(nsfn, deCaption))


tex_file <- file.path(opt$outdir,
                      paste(c("characterise.degenes",file_suffix,"tex"),collapse="."))

writeTex(tex_file, tex)
