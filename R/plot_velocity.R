stopifnot(
  require(ggplot2),
  require(reshape2),
  require(velocyto.R),
  require(Matrix),
  require(optparse),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--rdimstable"),
                default="none",
                help="A table containing the reduced coordinates and phenotype information"),
    make_option(c("--rdim1"),
                default="tSNE_1",
                help="The name of the column for reduced dimension one"),
    make_option(c("--rdim2"),
                default="tSNE_2",
                help="The name of the column for reduced dimension two" ),
    make_option(c("--matrixdir"),
                default=NULL,
                help=paste("Directory with matrix.mtx, barcodes.tsv and genes.tsv.",
                           "The matrix should contain exon (ex_ row nameprefix)",
                           "and intron (in_ row name prefix) counts" )),
    make_option(c("--minmaxclustavemat"), type="double",
                default=0.2,
                help="min.max.cluster.average for exon matrix filtering"),
    make_option(c("--minmaxclustavnmat"), type="double",
                default=0.05,
                help="min.max.cluster.average for intron matrix filtering"),
    make_option(c("--deltat"), type="double",
                default=1,
                help="deltaT parameter"),
    make_option(c("--kcells"), type="integer",
                default=15,
                help="number of cells for knn graph used for velocity estimates"),
    make_option(c("--fitquantile"), type="double",
                default=0.02,
                help="fit.quantile parameter"),
    make_option(c("--neighbourhoodsize"), type="integer",
                default=100,
                help="neighbourhood size"),
    make_option(c("--velocityscale"), type="character",
                default="log",
                help="velocity scale to use (log, sqrt, rank or linear)"),
    make_option(c("--arrowscale"), type="double",
                default=3,
                help="arrow scale multiplier"),
    make_option(c("--arrowlwd"), type="double",
                default=1.25,
                help="arrow line weight"),
    make_option(c("--gridflow"), type="logical",
                default=TRUE,
                help="show grid velocity summary"),
    make_option(c("--mingridcellmass"), type="double",
                default=0.5,
                help="min.grid.cell.mass parameter"),
    make_option(c("--gridn"), type="integer",
                default=40,
                help="number of grid points along each axis"),
    make_option(c("--cellalpha"), type="double",
                default=0.8,
                help="alpha value for cell border"),
    make_option(c("--cellborderalpha"), type="double",
                default=0,
                help="alpha value for cell border"),
    make_option(c("--showaxes"), type="logical",
                default=TRUE,
                help="show plot axes"),
    make_option(c("--plotdirvar"),
                default="velocityDir",
                help="latex var holding location of plots"),
    make_option(c("--plotcex"), type="double",
                default=1,
                help="overall cex for the plot"),
    make_option(c("--ncores"), type="integer",
                default=12,
                help="number of cores for velocity estimation"),
    make_option(c("--outdir"),
                default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

## read in the reduced dimension coordinates
rdims <- read.table(opt$rdimstable,
                    header=T, sep="\t",as.is=T)

print(dim(rdims))
print(head(rdims))

## read in the exon and intron count matrix
dat <- readMM(file.path(opt$matrixdir, "matrix.mtx"))
rownames(dat) <- read.table(file.path(opt$matrixdir,"genes.tsv"))$V1
colnames(dat) <- read.table(file.path(opt$matrixdir,"barcodes.tsv"))$V1

print(dim(dat))

## get the exon and intron counts for the relevant cells.
emat <- dat[grepl("ex_",rownames(dat)),
            colnames(dat) %in% rdims$barcode,drop=F]
nmat <- dat[grepl("in_",rownames(dat)),
            colnames(dat) %in% rdims$barcode,drop=F]

## trim off the exon and intron identifiers.
rownames(emat) <- substr(rownames(emat), 4, nchar(rownames(emat)))
rownames(nmat) <- substr(rownames(nmat), 4, nchar(rownames(nmat)))

## get the clusters
cluster.label <- factor(as.character(rdims$cluster))
names(cluster.label) <- rdims$barcode
nclust <- length(unique(rdims$cluster))

## set the cell colors according to cluster
cmap <- gg_color_hue(nclust)
names(cmap) <- c(0:(nclust-1))
cell.colors <- cmap[as.character(rdims$cluster)]
names(cell.colors) <- rdims$barcode

## set up the embedding
emb <- rdims[,c(opt$rdim1,opt$rdim2)]
rownames(emb) <- rdims$barcode

print(dim(emat))
print(dim(nmat))
print(head(emat))
print(head(nmat))
print(head(cluster.label))
print(table(cluster.label))

print(length(intersect(colnames(emat),names(cluster.label))))
print(length(intersect(colnames(nmat),names(cluster.label))))

## Filter genes by cluster expression
message("filtering emat")
emat <- filter.genes.by.cluster.expression(emat,
                                           cluster.label,
                                           min.max.cluster.average = opt$minmaxclustavemat)

message("filtering nmat")
nmat <- filter.genes.by.cluster.expression(nmat,
                                           cluster.label,
                                           min.max.cluster.average = opt$minmaxclustavnmat)

message("Number of genes with exon and intron counts:")
print(length(intersect(rownames(emat),rownames(emat))))

fit.quantile <- opt$fitquantile

# using non recommended default distance
rvel.cd <- gene.relative.velocity.estimates(emat,
                                            nmat,
                                            deltaT=opt$deltat,
                                            kCells=opt$kcells,
                                            fit.quantile=fit.quantile,
                                            n.cores=opt$ncores)

drawVelocityPlot <- function()
    {
        show.velocity.on.embedding.cor(emb=as.matrix(emb),
                                       vel=rvel.cd,
                                       n=opt$neighbourhoodsize,
                                       scale=opt$velocityscale,
                                       cell.colors=ac(cell.colors,alpha=opt$cellalpha),
                                       cex=opt$plotcex,
                                       arrow.scale=opt$arrowscale,
                                       arrow.lwd=opt$arrowlwd,
                                       show.grid.flow=opt$gridflow,
                                       min.grid.cell.mass=opt$mingridcellmass,
                                       grid.n=opt$gridn,
                                       do.par=F,
                                       cell.border.alpha = opt$cellborderalpha,
                                       axes=opt$showaxes)
    }

plotfilename="velocity.plot"

save_plots(file.path(opt$outdir, plotfilename),
           drawVelocityPlot,
           width=8,
           height=8)

texCaption <- paste("RNA velocity plot.")

tex <- paste(getSubsubsectionTex(texCaption),
             getFigureTex(plotfilename,texCaption,plot_dir_var=opt$plotdirvar),
             sep="\n")

print("Writing out latex snippet")

## write out latex snippet
tex_file <- file.path(opt$outdir,
                      "plot.velocity.tex")

writeTex(tex_file, tex)

message("Completed")
