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
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.Robj",
                help="A seurat object after PCA"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

outPrefix <- file.path(opt$outdir,"markers.summary")

cat("Running with options:\n")
print(opt)

s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)
s@ident <- cluster_ids

idents.all = sort(unique(cluster_ids))

if(!length(idents.all) > 1)
{
    stop("Only one cluster present")
}

gde.all = data.frame()

for (i in 1:length(idents.all)) {

    print(paste("Reading in cluster:",i))

    id <- idents.all[i]

    tableName = paste("markers.cluster",id,"txt","gz",
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

print(dim(filtered_markers))

## Annotate the filtered_markers with:
## (i) per-cluster expression levels
## (ii) per-cluster expression frequencies
## With lots of markers this is slow and memory intensive.
clusters <- as.numeric(unique(as.vector(cluster_ids)))
cluster_levels <- clusters[order(clusters)]
for(x in cluster_levels)
{
    print(paste("Adding expression information to the marker table for cluster:", x))
    fgenes <- filtered_markers$gene
    clust_cells <- names(cluster_ids[cluster_ids==x])
    other_cells <- names(cluster_ids[cluster_ids!=x])

    xmean <- apply(expm1(s@data[fgenes,clust_cells]),1,mean)
    omean <- apply(expm1(s@data[fgenes,other_cells]),1,mean)
    xfreq <- apply(s@data[fgenes,clust_cells],1,function(x) length(x[x>0])/length(x))

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
    xrows = filtered_markers$cluster==cluster
    tmp <- data.frame(filtered_markers[xrows])
    other <- cluster_levels[!cluster_levels==cluster]
    ocols <- paste(paste0("X",other),"exprs",sep="_")
    this <- paste(paste0("X",cluster),"exprs",sep="_")

    fcs <- log((tmp[,this]+1) / (tmp[,ocols]+1))

    ## for positive markers we store the min.
    min_fcs <- apply(fcs,1,min)

    ## for neg markers we store the max.
    max_fcs <- apply(fcs,1,max)

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
## differentially expressed marker genes as a txt file.

message("Saving marker.summary.table.txt")

marker_file <- paste(outPrefix,"table","txt","gz",
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

gp <- DoHeatmap(s, genes.use = top20$gene,
                slim.col.label = TRUE,
                remove.key = TRUE,
                cex.row=4,
                ## order.by.ident = TRUE, # depreciated
                title="Top 20 marker genes for each cluster")

save_ggplots(paste(outPrefix,"heatmap", sep="."),
           gp,
           width = 7,
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
            paste(outPrefix,"stats","txt",
                  sep="."),
            quote=F,sep="\t",row.names=F)
