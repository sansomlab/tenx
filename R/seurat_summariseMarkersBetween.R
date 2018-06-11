## Summarise the within-cluster differential expression analysis results
## for multiple clusters.
## run only specified comparison (in order for parallel execution)

# Libraries ----

stopifnot(
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(gplots),
  require(reshape2),
  require(xtable),
  require(ggplot2),
  require(openxlsx),
  require(optparse),
  require(tenxutils)
)

# Options ----

## deal with the options #mindiffpct, minpct, thresh.use
option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--testfactor"),default="factor",
                help="A column of @meta.data containing the conditions to be contrasted"),
    make_option(c("--a"), default="none",
                help="condition A"),
    make_option(c("--b"), default="none",
                help="condition B"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# set the run specs
run_specs <- paste(opt$numpcs,opt$resolution,opt$algorithm,opt$testuse,sep="_")

s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)

## cluster_ids
clusters <- sort(unique(as.vector(cluster_ids)))

acells <- s@meta.data[[opt$testfactor]] == opt$a
bcells <- s@meta.data[[opt$testfactor]] == opt$b

aName <- paste(opt$a,"mean",sep="_")
bName <- paste(opt$b,"mean",sep="_")

ntake <- 20

start <-T
begin <- T
for(cluster in clusters)
{
    as <- names(cluster_ids)[acells & cluster_ids==cluster]
    bs <- names(cluster_ids)[bcells & cluster_ids==cluster]
    no.a <- length(as)
    no.b <- length(bs)
    if(start)
    {
        stats <- data.frame(cluster=c(cluster),
                            no.a=c(no.a),
                            no.b=c(no.b),stringsAsFactors=F)
        start <- F
    } else {
        stats <- rbind(stats, c(as.numeric(cluster), no.a, no.b))
        }

    ## get the differentially expressed genes
    res_fn <- file.path(
        opt$outdir,
        paste("markers.between",opt$testfactor,
              "cluster",cluster,
              "txt","gz",
              sep="."))
    message("res_fn:")
    print(res_fn)

    ## Note that we expect valid cases where res_fn does not exist
    ## i.e. when marker genes were not identified.

    if(file.exists(res_fn))
       {

           de <- read.table(gzfile(res_fn), header=T, sep="\t", as.is=T)
           de <- de[order(de$p_va),]

           ameans <- rowMeans(s@scale.data[de$gene,as])
           bmeans <- rowMeans(s@scale.data[de$gene,bs])

           de[[aName]] <- ameans
           de[[bName]] <- bmeans

           detop <- de[1:min(nrow(de),ntake),]

           if(begin)
           {
               res <- de
               high <- detop
               begin <- F
           }
           else {
               res <- rbind(res, de)
               high <- rbind(detop,high)
           }
       }
}

res <- res[,c("cluster","gene","p.adj","p_val","avg_logFC",
              "pct.1","pct.2",aName,bName,"gene_id")]

out_fn <- file.path(
    opt$outdir,
    paste("markers.between",opt$testfactor,
          "summary.table.txt.gz",sep=".")
)

## write out the full table of differentially expressed genes.
## in txt format
message("Saving markers to txt file")
write.table(res, gzfile(out_fn), quote=F,
            row.names=F, col.names=T,
            sep="\t")

## write out a p-value filtered table (max 10% FDR) in
## in excel format
message("Saving markers in excel format")

xres <- res[res$p.adj < 0.1,]

wb <- createWorkbook()

addWorksheet(wb,"markers_between")
setColWidths(wb,"markers_between",cols=1:ncol(xres),widths=10)
hs <- createStyle(textDecoration = "BOLD")
writeData(wb, "markers_between", tidyNumbers(xres),
          withFilter = T, headerStyle=hs)
saveWorkbook(wb, file=gsub(".txt.gz",".xlsx",out_fn),
             overwrite=T)



##
rownames(stats) <- stats$cluster
stats$nde <- 0
for(cluster in clusters)
{
    nde <- nrow(res[res$p.adj < 0.1 & res$cluster==cluster,])
    stats[cluster,"nde"] <- nde

}

colnames(stats) <- c("cluster",
                     paste0("ncells_",opt$a),
                     paste0("ncells_",opt$b),
                     "n_de_genes")

tex_fn <- file.path(
    opt$outdir,
    paste("markers.between.summary.tex",sep=".")
)

xtab <- stats
xtab$cluster <- NULL
print(xtable(xtab, digits=0,
             caption="summary of clusters by condition"),
      file=tex_fn)

pf <- melt(stats[,1:3],id.vars="cluster")

gp <- ggplot(pf,aes(cluster,value,fill=variable))
gp <- gp + geom_bar(stat="identity",position="dodge")
gp <- gp + ylab("number of cells")
gp <- gp + scale_fill_manual(values=c("seagreen4","bisque2"))

save_ggplots(gsub(".tex",".cellnumbers",tex_fn),
             gp,
             width=7,
             height=7)

gp <- ggplot(stats,aes(cluster,n_de_genes))
gp <- gp + geom_bar(stat="identity")
gp <- gp + ylab("number of de genes")

save_ggplots(gsub(".tex",".nde",tex_fn),
             gp,
             width=5,
             height=5)

## Make a heatmap
diffMat <- dcast(res, gene~cluster, value.var="avg_logFC")
diffMat[is.na(diffMat)] <- 0

rownames(diffMat) <- diffMat$gene
diffMat$gene <- NULL


m <- as.matrix(diffMat[unique(high$gene),])

ramp_colors <- c("blue","cyan","black","yellow","red")

nbreaks=100 # number of graduations
rm <- range(m)
rm <- c(-1.2,1.2)
breaks=seq(rm[1],rm[2],diff(rm)/nbreaks) #range (-2 -> 2) and breakpoint size = range_size / no. graduations
colors=colorRampPalette(ramp_colors)(nbreaks) # get the color pallete


plotfn <- gsub(".tex",".heatmap",tex_fn)

plot_fn <- function()
    {
        par(cex.main=0.7)
        heatmap.2(m,
                  col=colors,
                  breaks=breaks,
                  scale="none",
                  Colv=F,
                  Rowv=F,
                  mar=c(4,7),
                  trace="none",
                  key.title = "",
                  key.par=list(mar=c(2,0.3,2,0.3)),
                  density.info=c("none"),
                  lwid = c(1,5),
                  lhei = c(1,8),
                  key.xlab = "avg_logFC",
                  key.ylab = "",
                  xlab="cluster",
                  cexRow = 0.4,
                  cexCol = 1.3
                  )
    }

save_plots(plotfn,
           plot_fn=plot_fn,
           width=6,
           height=8)
