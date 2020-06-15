## Summarise the marker genes across a set of clusters

# Libraries ----

stopifnot(
  require(Matrix),
  require(reshape2),
  require(data.table),
  require(optparse)
)

# Options ----

option_list <- list(
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

outPrefix <- file.path(opt$outdir,"markers.summary")
if(!is.null(opt$subgroup)) { opt$subgroup <- strsplit(opt$subgroup,",")[[1]]}
cat("Running with options:\n")
print(opt)
cluster_ids <- readRDS(opt$clusterids)
idents.all = sort(unique(cluster_ids))

if(!length(idents.all) > 1)
{
    stop("Only one cluster present")
}

## Construct a data frame containing the findMarkers results for all clusters
begin <- TRUE
for (i in 1:length(idents.all)) {

    print(paste("Reading in cluster:",i))

    id <- idents.all[i]

    tableName = paste("cluster.stats",id,"txt","gz",
                      sep=".")

    markerFile <- file.path(opt$outdir,tableName)
    message("markerFile: ", markerFile)

    if(file.exists(markerFile))
        {
            stats <- read.table(gzfile(markerFile),
                                header=T,
                                row.names=1,
                                as.is=T,
                                sep="\t")

            colnames(stats) <- paste0("x", id, "_", colnames(stats)
                                      )
            if(begin) {
                stats.all <- stats
                begin <- FALSE

            } else {
                if(rownames(stats) == rownames(stats.all))
                   {
                       stats.all = cbind(stats.all, stats)
                   } else {
                       stop("clusters stats have different rownames...")
                   }

            }

        } else {
            print(paste("No marker file found for cluster",i,sep=" "))
        }

}

stats.all$gene <- rownames(stats.all)

out_path = file.path(opt$outdir, "cluster.stats.summary.table.txt.gz")
write.table(stats.all, gzfile(out_path),
            sep="\t", col.names=T, row.names=T, quote=F)


message("completed")
