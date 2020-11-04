## Summarise the marker genes across a set of clusters

# Libraries ----

stopifnot(
  require(Seurat),
  require(dplyr),
  require(Matrix),
  require(reshape2),
  require(data.table),
  require(openxlsx),
  require(optparse),
  require(ComplexHeatmap),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--seuratassay"), default="RNA",
                help="The seurat assay to use"),
    make_option(c("--statstable"), default="none",
                help="The table of per-cluster statistics"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--subgroup"), default=NULL,
                help="Optional. Name of a column in metadata to add annotation for in the heatmap"),
    make_option(c("--pdf"), default = FALSE,
                help="Create a pdf version of the top marker heatmap"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

outPrefix <- file.path(opt$outdir,"markers.summary")
if(!is.null(opt$subgroup)) { opt$subgroup <- strsplit(opt$subgroup,",")[[1]]}
cat("Running with options:\n")
print(opt)

s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
Idents(s) <- cluster_ids

message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("seurat_summariseMarkers.R running with default assay: ", DefaultAssay(s))

idents.all = sort(unique(cluster_ids))

if(!length(idents.all) > 1)
{
    stop("Only one cluster present")
}

## Construct a data frame containing the findMarkers results for all clusters
gde.all = data.frame()

for (i in 1:length(idents.all)) {

    print(paste("Reading in cluster:",i))

    id <- idents.all[i]

    tableName = paste("markers.cluster",id,"tsv","gz",
                      sep=".")

    markerFile <- file.path(opt$outdir,tableName)
    message("markerFile: ", markerFile)

    if(file.exists(markerFile))
        {
            gde <- read.table(gzfile(markerFile),
                              header=T,
                              as.is=T,
                              sep="\t")

            print(paste("Length of cluster marker matrix:",dim(gde)[1]))
            gde.all = rbind(gde.all, gde)
        } else {
            print(paste("No marker file found for cluster",i,sep=" "))
        }

}

markers <- gde.all

print("Marker aggregation complete")
print(dim(markers))

markers <- markers[,c("cluster","gene","gene_id",
                      "p_val","p.adj",
                      "avg_logFC","pct.1","pct.2",
                      "cluster_mean","other_mean")]

markers <- markers[order(markers$cluster, markers$p.adj),]

markers$index <- c(rownames(markers))

print(dim(markers))

print("Filtering by adjusted p-value")
filtered_markers <- data.table(markers[markers$p.adj < 0.1,])
message("Filtered markers are:")
print(dim(filtered_markers))

## Annotate the filtered_markers with:
## (i) per-cluster expression levels
## (ii) per-cluster expression frequencies
## With lots of markers this is slow and memory intensive.
stats <- read.table(gzfile(opt$statstable), sep="\t", header=T,
                    as.is=T)
rownames(stats) <- stats$gene

print(head(stats))

clusters <- as.numeric(unique(as.vector(cluster_ids)))
cluster_levels <- clusters[order(clusters)]
for(x in cluster_levels)
{
    # skip the catch-all 911 cluster.
    if(as.character(x) == "911") { next }

    print(paste("Adding expression information to the marker table for cluster:", x))
    fgenes <- filtered_markers$gene
    if(!all(fgenes %in% rownames(stats)))
    {
        stop("cluster stats table does not contain all of the genes")
    }

    clust_cells <- names(cluster_ids[cluster_ids==x])
    other_cells <- names(cluster_ids[cluster_ids!=x])

    # xmean <- apply(expm1(GetAssayData(object = s, slot="data")[fgenes,clust_cells]),1,mean)
    # omean <- apply(expm1(GetAssayData(object = s, slot="data")[fgenes,other_cells]),1,mean)
    # xfreq <- apply(GetAssayData(object = s, slot="data")[fgenes,clust_cells],1,function(x) length(x[x>0])/length(x))
    ## xmean <- Matrix::rowMeans(expm1(GetAssayData(object = s, slot="data")[fgenes,clust_cells]))
    ## omean <- Matrix::rowMeans(expm1(GetAssayData(object = s, slot="data")[fgenes,other_cells]))
    ## xfreq <- Matrix::rowSums(GetAssayData(object = s, slot="data")[fgenes,clust_cells]>0)/length(clust_cells)

    xmean <- expm1(stats[fgenes, paste0("x",x,"_cluster_mean")])
    omean <- expm1(stats[fgenes, paste0("x",x,"_other_mean")])
    xfreq <- stats[fgenes,paste0("x",x,"_cluster_freq")]

    filtered_markers[[paste0(x,"_exprs")]] <- xmean
    filtered_markers[[paste0(x,"_freq")]] <- xfreq

    ## we add 0.1 to the mean counts to moderate the fold change
    ## - a count of 1 in 33% of cells vs 0 = 4 fold change
    ## - a count of 1 in 33% of cells vs 1 in 10% = 2 fold change

    xrows <- filtered_markers$cluster==x
    filtered_markers$mLog2FC[xrows] <- log2((xmean[xrows]+0.1)/(omean[xrows]+0.1))
    filtered_markers$exprs[xrows] <- xmean[xrows]
    filtered_markers$freq[xrows] <- xfreq[xrows]
}

## compute the minimum fold change to another cluster for the filtered markers
filtered_markers$min_logFC <- NA
cluster_levels <- unique(filtered_markers$cluster)

for(cluster in cluster_levels)
{
    ## skip the catch-all 911 cluster.
    if(as.character(cluster) == "911") { next }

    message("computing min fold changes for cluster:", cluster)
    xrows = filtered_markers$cluster==cluster
    tmp <- data.frame(filtered_markers[xrows])
    other <- cluster_levels[!cluster_levels==cluster]
    ocols <- paste(paste0("X",other),"exprs",sep="_")
    this <- paste(paste0("X",cluster),"exprs",sep="_")

    this_exprs <- tmp[,this]+1
    other_exprs <- tmp[,ocols]+1
    fcs <- log(this_exprs / other_exprs )

    if(length(ocols) > 1)
    {
        ## for positive markers we store the min.
        min_fcs <- apply(fcs,1,min)

        ## for neg markers we store the max.
        max_fcs <- apply(fcs,1,max)
    } else {
        # there is only 1 other column so
        min_fcs <- max_fcs <- fcs
    }

    filtered_markers$min_logFC[xrows] <- min_fcs
    filtered_markers$max_logFC[xrows] <- max_fcs
}

markers$min_logFC <- NA
markers$max_logFC <- NA

## rownames(markers) <- paste(markers$cluster,markers$gene_id,sep="-")
## rownames(filtered_markers) <- paste(filtered_markers$cluster,filtered_markers$gene_id,sep="-")
take_cols <- c("min_logFC","max_logFC")
markers[filtered_markers$index, take_cols] <- filtered_markers[,take_cols,with=F]

## write out the full (unannotated) table of significantly
## differentially expressed marker genes as a.tsv file.

message("Saving marker.summary.table.tsv")

marker_file <- paste(outPrefix,"table","tsv","gz",
                     sep=".")

write.table(markers,
            gzfile(marker_file),
            quote=F,sep="\t",row.names=F)


## write out a shorter table of significantly differentially
## expressed marker genes in excel format.
## TODO allow specification of adjusted p-value threshold.

message("Saving marker.summary.table.xlsx")
wb <- createWorkbook()


addWorksheet(wb,"filtered_markers")
setColWidths(wb,"filtered_markers",cols=1:ncol(filtered_markers),widths=10)
hs <- createStyle(textDecoration = "BOLD")
writeData(wb, "filtered_markers", tidyNumbers(filtered_markers),
          withFilter = T, headerStyle=hs)
saveWorkbook(wb, file=paste(outPrefix,"table","xlsx",
                            sep="."),
             overwrite=T)

message("Making a heatmap of the top marker genes from each cluster")

## make a heatmap of the top DE genes.
filtered_markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20

if(!is.null(opt$subgroup))
{
    if(!opt$subgroup %in% colnames(s@meta.data))
    {
        opt$subgroup <- NULL
    }
}

## y <- x[x$cluster==0,]

## only draw the plot if scale.data is populated.

if(max(dim(GetAssayData(s, slot="scale.data"))) > 0 )
{
    mch <- markerComplexHeatmap(s,
                                marker_table=filtered_markers,
                                n_markers=20,
                                cells_use=NULL,
                                row_names_gp=11,
                                sub_group=opt$subgroup)

    drawHeatmap <- function()
    {
        draw(mch)
    }

} else {

    drawHeatmap <- function()
    {
    plot.new()
    text(0.5,0.5,"scale.data slot not present")
    }
}



save_plots(paste(outPrefix,"heatmap", sep="."),
           plot_fn=drawHeatmap,
           width = 7,
           to_pdf = opt$pdf,
           height = 9)



## summarise the number of marker genes identified for each cluster
summary <- c()
for(id in idents.all)
{
    ncells = length(cluster_ids[cluster_ids==id])
    npos = length(filtered_markers$p.adj[filtered_markers$cluster==id & filtered_markers$avg_logFC > 0] )
    nneg = length(filtered_markers$p.adj[filtered_markers$cluster==id & filtered_markers$avg_logFC < 0] )
    ntotal = npos + nneg
    summary <- c(summary,c(id, ncells, npos, nneg, ntotal))
}

sumdf <- data.frame(matrix(summary,ncol=5,byrow=T))
colnames(sumdf) <-c("cluster","ncells","n_pos_markers","n_neg_markers","n_markers")

print("Writing out some summary statistics")
write.table(sumdf,
            paste(outPrefix,"stats","tsv",
                  sep="."),
            quote=F,sep="\t",row.names=F)


message("seurat_summariseMarkers.R final default assay: ", DefaultAssay(s))

message("completed")
