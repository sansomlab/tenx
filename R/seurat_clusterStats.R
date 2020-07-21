## Compute statistics for a given cluster
## Note that p-values must be adjusted in down-stream analysis
## when used for the analysis of multiple clusters.
## run only specified comparison (in order for parallel execution)

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(SeuratDisk),
  require(future),
  require(dplyr),
  require(Matrix),
  require(tenxutils),
  require(reshape2)
)


# Options ----

option_list <- list(
    make_option(c("--seuratobject"), default="begin.rds",
                help="A seurat object after PCA"),
    make_option(c("--seuratassay"), default="RNA",
                help="the assay to set as the default"),
    make_option(c("--clusterids"), default="none",
                help="A list object containing the cluster identities"),
    make_option(c("--cluster"), default=1,
                help="The identity of the cluster to test"),
    make_option(c("--testfactor"),default="none",
                help="A column of @meta.data containing the conditions to be contrasted"),
    make_option(c("--a"), default="none",
                help="condition A"),
    make_option(c("--b"), default="none",
                help="condition B"),
    make_option(c("--conservedfactor"), default="none",
                help="A column of @meta.data containing a grouping factor across which markers should be conserved"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir"),
        make_option(
        c("--numcores"),
        type="integer",
        default=1,
        help="Number of cores to be used for the Jackstraw analysis"
        )

    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


s <- loadSeurat(path=opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)

xx <- names(cluster_ids)
yy <- colnames(x=s)
if(!all.equal(xx[order(xx)], yy[order(yy)])) {
    stop("cluster_ids and seurat object have different cells")
}

cluster_ids <- cluster_ids[colnames(x = s)]

if(!identical(colnames(x = s),names(cluster_ids)))
{   stop("Cluster cell names do not match Seurat object cell names")
}

if(!identical(names(cluster_ids), rownames(s[[]])))
{
    # probaby not necessary.
    stop("cluster_ids and metadata rownames do not match")
}


# Set the identities of the cells groups to be tested
id <- opt$cluster #idents.all[i]
ident <- rep("ignore",length(cluster_ids))
names(ident) <- names(cluster_ids)
# Set the grouping factor across which markers should be conserved
ident.conserved <- factor(rep("all",length(cluster_ids)))
names(ident.conserved) <- names(cluster_ids)

if(opt$testfactor=="none")
{
    ident[cluster_ids == id] <- "a"
    ident[cluster_ids != id] <- "b"
    if(length(ident[ident == "ignore"])>0)
    { stop("Problem assigning cells to clusters") }
} else {
     print(head(s[[opt$testfactor]]))
    if(!opt$a %in% s[[]][,opt$testfactor] | !opt$b %in% s[[]][,opt$testfactor])
    {
        stop("between_a and/or between_b not present in the opt$testfactor metadata colum")
    }
    ident[cluster_ids==id & s[[]][,opt$testfactor] == opt$a] <- "a"
    ident[cluster_ids==id & s[[]][,opt$testfactor] == opt$b] <- "b"
}

if(opt$conservedfactor != "none"){
    # stopifnot testfactor also supplied
    ident.conserved <- factor(s[[]][names(ident), opt$conservedfactor])
    if (nlevels(ident.conserved) < 2){
        stop("Conserved factor has fewer than 2 levels")
    }
}

message("Setting default assay to: ", opt$seuratassay)
DefaultAssay(s) <- opt$seuratassay

message("seurat_FindMarkers.R running with default assay: ", DefaultAssay(s))

stats.list <- list()

# one level ('ignore') if no conservation is required
# as many levels as
for (conserved.level in levels(ident.conserved)){
    message("conserved.level: ", conserved.level)
    # Restore cluster identity for all cells
    ident.use <- ident
    # Ignore cells in other levels of conservation factor
    ident.use[ident.conserved != conserved.level] <- "ignore"
    Idents(s) <- factor(ident.use)
    ident.use = Idents(s)
    message("Identities:")
    print(table(Idents(s)))

    # Original FindMarker pipeline routine
    ## check fold change threshold
    nclust_check <- F
    ncells_check <- F
    other_checks <- F

    if(length(unique(Idents(s)[Idents(s)!="ignore"])) > 1)
    {
        nclust_check = TRUE

        cluster_cells <- colnames(x = s)[Idents(s) == "a"]
        other_cells <- colnames(x = s)[Idents(s) == "b"]

        ## compute percentages and difference
        genes <- rownames(x = s)

	cluster_pct <- rowSums(GetAssayData(object = s, slot="data")[genes,cluster_cells, drop=F]>0)/length(cluster_cells)
	other_pct <- rowSums(GetAssayData(object = s, slot="data")[genes,other_cells, drop=F]>0)/length(other_cells)

        pcts <- cbind(cluster_pct,other_pct)

	message("calc mean expression in cluster...")
	cluster_mean <- FastExpMeanChunked(GetAssayData(object = s, slot="data")[genes,cluster_cells],2000)
	message("calc mean expression in other cells...")
	other_mean <- FastExpMeanChunked(GetAssayData(object = s, slot="data")[genes,other_cells],2000)

	message("saving stats")

        ## store these stats so that they can added to the results table
        ## e.g. for later investigation of threshold effects...

        cluster_stats <- data.frame(row.names=rownames(x = s),
                                   cluster_mean=signif(cluster_mean,4),
                                   other_mean=signif(other_mean,4),
                                   cluster_pct=signif(cluster_pct,4),
                                   other_pct=signif(other_pct,4))


        if (opt$conservedfactor != "none"){
            stats.list[[conserved.level]] <- cluster_stats
        }

        print(head(cluster_stats))
    } else {
        print("Only one cluster/condition!")
    }
}



if (length(stats.list) > 1){
    message("Multiple sets of markers were detected. Identifying conserved markers ...")

    conserved.table <- data.frame(
        row.names=rownames(x = s),
        cluster = rep(opt$cluster, length(rownames(x = s)))
    )

    message("conserved table initialised")

    for(stat_name in colnames(cluster_stats))
        {

            ## average percentage of cells in group 1 ----
            tmp.table <- do.call(
                "cbind",
                lapply(stats.list, function(x){
                    x[, stat_name]
                })
            )
            conserved.table[[stat_name]] <- rowMeans(as.matrix(tmp.table))
        }

    print(head(conserved.table))

    cluster_stats <- conserved.table

}


out_path <- file.path(
    opt$outdir,
    paste("cluster.stats",opt$cluster,"tsv","gz",sep="."))

print(paste("Saving markers to:", out_path))

write.table(cluster_stats,
            gzfile(out_path),
            quote=F,sep="\t",row.names=T)


message("seurat_ClusterStats.R final default assay: ", DefaultAssay(s))

message("completed")
