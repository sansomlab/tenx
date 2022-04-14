## Identify marker genes for a given cluster.
## Note that p-values must be adjusted in down-stream analysis
## when used for the analysis of multiple clusters.
## run only specified comparison (in order for parallel execution)

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
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
    make_option(c("--testuse"), default="wilcox",
                help="test to use for finding cluster markers"),
    make_option(c("--minpct"), default=0.1,
                help="minimum fraction of cells expressing gene"),
    make_option(c("--mindiffpct"), default=-Inf,
                help="minimum fraction difference between cluster and other cells. Setting not recommended."),
    make_option(c("--threshuse"), default=0.25,
                help="testing limited to genes with this (log scale) difference in mean expression level."),
    make_option(c("--mincells"), default=3,
                help="minimum number of cells required (applies to cluster and to the other cells)"),
    make_option(c("--maxcellsperident"), default=Inf,
                help="max.cell.per.ident"),
    make_option(c("--testfactor"),default="none",
                help="A column of @meta.data containing the conditions to be contrasted"),
    make_option(c("--a"), default="none",
                help="condition A"),
    make_option(c("--b"), default="none",
                help="condition B"),
    make_option(c("--conservedfactor"), default="none",
                help="A column of @meta.data containing a grouping factor across which markers should be conserved"),
    make_option(c("--conservedpadj"), default=0.05,
                help="adjusted P-value cutoff required in all individual tests for conserved markers"),
    make_option(c("--annotation"), default="not_set",
                help=paste('A gzipped file with "ensembl_id" and "gene_name" columns.',
                           'Used if s@misc$gene_id is missing')),
    make_option(c("--universe"), default=F,
                help="write out the gene universe (e.g. for GSEA)"),
    make_option(c("--testname"), default="default",
                help="the name for this test"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
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

plan("multiprocess",
     workers = opt$numcores)

options(future.globals.maxSize=100000000000000000)

s <- readRDS(opt$seuratobject)
cluster_ids <- readRDS(opt$clusterids)[Cells(s)]

have_gene_ids <- "gene_id" %in% colnames(s@misc)

if(!have_gene_ids)
{
    cat("Missing ensembl gene_id column in @misc, reading in annotation.\n")

    ann <- read.table(gzfile(opt$annotation), header=T, sep="\t", as.is=T)

    ## subset to columns of interest
    ann <- unique(ann[,c("ensembl_id", "gene_name")])

    ## the gene names have been made unique so we want only 1:1 mappings.
    gn_tab <- table(ann$gene_name)

    unique_gn <- names(gn_tab)[gn_tab == 1]
    ann <- ann[ann$gene_name %in% unique_gn,]

    rownames(ann) <- ann$gene_name
}

xx <- names(cluster_ids)
yy <- colnames(x=s)
if(!isTRUE(all.equal(xx[order(xx)], yy[order(yy)]))) {
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

markers.conserved.list <- list()
background.conserved.list <- list()

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

        if(length(cluster_cells) >= opt$mincells & length(other_cells) >= opt$mincells)
        {
            ncells_check <- TRUE
        } else {
            print("Either the cluster or the other cells have fewer than the minimum number")
        }

        ## compute percentages and difference
        genes <- rownames(x = s)
        #cluster_pct <- apply(GetAssayData(object = s)[genes,cluster_cells, drop=F], 1, function(x) round(sum(x>0)/length(x), digits=3))
        #other_pct <- apply(GetAssayData(object = s)[,other_cells, drop=F], 1, function(x) round(sum(x>0)/length(x),digits=3))
	cluster_pct <- round(rowSums(GetAssayData(object = s, slot="data")[genes,cluster_cells, drop=F]>0)/length(cluster_cells),digits=3)
	other_pct <- round(rowSums(GetAssayData(object = s, slot="data")[genes,other_cells, drop=F]>0)/length(other_cells),digits=3)

        pcts <- cbind(cluster_pct,other_pct)
        rm(cluster_pct)
        rm(other_pct)

        max_pct <- apply(pcts,1,max)
        min_pct <- apply(pcts,1,min)
        diff_pct <- max_pct - min_pct

        rm(min_pct)
        rm(pcts)

        ## compute mean expression levels and difference
        #cluster_mean <- apply(GetAssayData(object = s)[,cluster_cells, drop=F], 1, FUN = ExpMean)
        #other_mean <- apply(GetAssayData(object = s)[,other_cells, drop=F], 1, FUN = ExpMean)
        #diff_mean <- abs(cluster_mean - other_mean)

	message("calc mean expression in cluster...")
	cluster_mean <- FastExpMeanChunked(GetAssayData(object = s, slot="data")[genes,cluster_cells],2000)
	message("calc mean expression in other cells...")
	other_mean <- FastExpMeanChunked(GetAssayData(object = s, slot="data")[genes,other_cells],2000)
	diff_mean <- abs(cluster_mean - other_mean)
	message("saving stats")

        ## store these stats so that they can added to the results table
        ## e.g. for later investigation of threshold effects...

        filter_stats <- data.frame(row.names=rownames(x = s),
                                   cluster_mean=signif(cluster_mean,4),
                                   other_mean=signif(other_mean,4))

        rm(cluster_mean)
        rm(other_mean)

        ## make an identity vector of genes satisfying the given criteria
        take <- max_pct > opt$minpct & diff_pct > opt$mindiffpct & diff_mean > opt$threshuse
        rm(diff_pct)
        rm(diff_mean)

        ## deliberately ignore diff_mean when selecting the background genes
        background_take <- max_pct > opt$minpct
        rm(max_pct)

        genes.use <- rownames(x = s)[take]
        background.use <- rownames(x = s)[background_take]

        if(length(genes.use > 0))
        {
            print(paste("Number of genes to be tested:",length(genes.use)))
            other_checks = TRUE
        } else {
            print("Other checks failed: no genes satisfied opt$minpct & opt$mindiffpct & opt$theshuse")
        }

    } else {
        print("Only one cluster/condition!")
    }

    if(nclust_check & ncells_check & other_checks)
    {

        ## Marker identification using customised FindMarkers routine
        ## handle ROC special case
        ## otherwise all p values kept for adjusting.
        if(opt$testuse=="ROC")
        {
            stop("ROC not supported")
        }

        ident.use = Idents(s)
        ##idents.all = sort(unique(Idents(s)))

        print(table(ident.use)) ##idents.all)
        #genes.de = list()

        min.cells = 3

        id <- opt$cluster #idents.all[i]

        markers <- FindMarkers(s,
                               ident.1 = "a",
                               ident.2 = "b",

                               logfc.threshold = opt$threshuse,
                               test.use = opt$testuse,
                               min.pct = opt$minpct,
                               min.diff.pct = opt$mindiffpct,

                               min.cells.feature = min.cells,
                               min.cells.group = min.cells,

                               max.cells.per.ident=opt$maxcellsperident,
                               base = exp(1))

        print(head(markers))
        ## keep everything, adjust later
        return.thresh = 1

        if (nrow(markers) > 0) {
            markers = markers[order(markers$p_val, -markers$avg_logFC), ]
            markers = subset(markers, p_val < return.thresh)

            if (nrow(markers) > 0){
                markers$cluster = opt$cluster
            }
            markers$gene = rownames(markers)
        }

        print("FindMarkers complete")

        message("correcting p-values")
        ## Add a BH corrected p-value
        ## correction is deliberately applied separately within cluster
        markers$p.adj <- p.adjust(markers$p_val, method="BH")

        message("selecting columns of interest")
        markers <- markers[,c("cluster","gene","p.adj","p_val","avg_logFC","pct.1","pct.2")]
        markers <- markers[order(markers$cluster, markers$p_val),]

        print(head(markers))
        print(dim(markers))

        ## add the ensembl gene_ids...
        if(have_gene_ids)
        {
            print("Adding enesmbl gene_ids from @misc")
            markers$gene_id <- s@misc[markers$gene,"gene_id"]
        } else {
            print("Adding ensembl gene_ids from annotation")
            markers$gene_id <- ann[markers$gene, "ensembl_id"]
        }


        message("adding the filter stats")
        ## add the filter stats for each gene
        markers <- cbind(markers, filter_stats[rownames(markers),])

        if (opt$conservedfactor != "none"){
            markers.conserved.list[[conserved.level]] <- subset(markers, p.adj < opt$conservedpadj)
        }

        ## write out the full results
        ## file name should be markers.cluster.x.tsv

        if(opt$testfactor=="none"){
            if (opt$conservedfactor == "none"){
                prefix = "markers.cluster"
            } else {
                prefix = sprintf("markers.%s.cluster", conserved.level)
            }
        }
        else {
            if (opt$conservedfactor == "none"){
                prefix = sprintf("markers.between.%s.cluster", opt$testfactor)
            } else {
                prefix = sprintf("markers.%s.between.%s.cluster", conserved.level, opt$testfactor)
            }
        }

        out_path <- file.path(
            opt$outdir,
            paste(prefix,opt$cluster,"tsv","gz",sep="."))

        print(paste("Saving markers to:", out_path))

        write.table(markers,
                    gzfile(out_path),
                    quote=F,sep="\t",row.names=F)

        rm(markers)
        ## write out the ensembl gene_ids of the "universe" for downstream
        ## geneset analysis
        ##
        ## * these are taken to be the genes passing min.pct
        ##   _but not_ min.diff.pct or thresh.use

        if (opt$conservedfactor != "none"){
            message("adding background list for level:")
            print(length(background.use))
            background.conserved.list[[conserved.level]] <- background.use
        }

        if(have_gene_ids)
        {
            print("Writing universe gene ids using @misc$gene_id")
            universe_gene_ids = s@misc[background.use,"gene_id"]
        } else {
            print("Writing universe gene ids using annotation file")
            # avoid writing out NA values.
            background.take = background.use[background.use %in% rownames(ann)]
            universe_gene_ids = ann[background.take,"ensembl_id"]
        }

        universe <- data.frame(gene_id = universe_gene_ids)

        universe_path <- file.path(
            opt$outdir,
            paste(prefix,opt$cluster,"universe","tsv","gz",sep="."))

        write.table(universe,
                    gzfile(universe_path),
                    quote=FALSE, sep="\t", row.names=FALSE)

    } else {
        print("Skipping find marker routine, one or more checks failed to pass")
    }

}

if (length(markers.conserved.list) > 1){
    message("Multiple sets of markers were detected. Identifying conserved markers ...")
    message("Individual count of markers:")
    print(sapply(markers.conserved.list, "nrow"))
    ## Debug: RDS lists of marker table for individual conservation levels
    # if(opt$testfactor=="none"){
    #     rdsFile <- file.path(opt$outdir, sprintf("markers.conserved.cluster.%s.list.rds", opt$cluster))
    # } else {
    #     rdsFile <- file.path(opt$outdir, sprintf("markers.between.conserved.cluster.%s.list.rds", opt$cluster))
    # }
    # saveRDS(markers.conserved.list, rdsFile)
    # Identify markers conserved across all conservation levels
    markers.sig <-
        table(unlist(lapply(markers.conserved.list, "rownames"), use.names = FALSE)) == length(markers.conserved.list)
    message("Markers significant in all levels:")
    table(markers.sig)
    markers.sig <- names(markers.sig)[markers.sig]

    markers.fc <- do.call("cbind", lapply(markers.conserved.list, function(x){x[markers.sig, "avg_logFC"]}))
    rownames(markers.fc) <- markers.sig
    colnames(markers.fc) <- names(markers.conserved.list)

    same.dir <- rowSums(abs(markers.fc)) == abs(rowSums(markers.fc))
    message("significant in the same direction:")
    print(table(same.dir))

    markers.conserved <- names(same.dir)[same.dir]
    message("markers conserved:")
    print(length(markers.conserved))

    # general information ----

    # avoid corner-case whereby no marker is conserved across all levels
    if (length(markers.conserved) > 0){

        conserved.table <- data.frame(
            cluster = opt$cluster,
            gene = markers.conserved.list[[1]][markers.conserved, "gene"],
            gene_id = markers.conserved.list[[1]][markers.conserved, "gene_id"],
            row.names = markers.conserved
        )

        message("conserved table initialised")

        # max p_val across levels ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "p_val"]
          })
        )
	print(str(tmp.table))
        rownames(tmp.table) <- markers.conserved
        conserved.table$p_val <- apply(
            as.matrix(tmp.table[markers.conserved, , drop=FALSE]), 1, max)

        message("max p_val added")

        # max p.adj across levels ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "p.adj"]
          })
        )
        rownames(tmp.table) <- markers.conserved
        conserved.table$p.adj <- apply(
            as.matrix(tmp.table[markers.conserved, , drop=FALSE]), 1, max)

        # average fold change across levels ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "avg_logFC"]
          })
        )
        rownames(tmp.table) <- markers.conserved
        tmp.table <- exp(tmp.table) - 1
        conserved.table$avg_logFC <- log(rowMeans(as.matrix(tmp.table[markers.conserved, , drop=FALSE])) + 1)

        message("avg_logFC added")

        # average percentage of cells in group 1 ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "pct.1"]
          })
        )
        rownames(tmp.table) <- markers.conserved
        conserved.table$pct.1 <- rowMeans(as.matrix(tmp.table[markers.conserved, , drop=FALSE]))

        # average percentage of cells in group 2 ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "pct.2"]
          })
        )
        rownames(tmp.table) <- markers.conserved
        conserved.table$pct.2 <- rowMeans(as.matrix(tmp.table[markers.conserved, , drop=FALSE]))

        message("pct of cells groups 1 and 2 added")

        # mean in level:cluster ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "cluster_mean"]
          })
        )
        rownames(tmp.table) <- markers.conserved
        tmp.table <- exp(tmp.table) - 1
        conserved.table$cluster_mean <- log(rowMeans(as.matrix(tmp.table[markers.conserved, , drop=FALSE])) + 1)

        # mean in level:other_clusters ----
        tmp.table <- do.call(
          "cbind",
          lapply(markers.conserved.list, function(x){
            x[markers.conserved, "other_mean"]
          })
        )
        rownames(tmp.table) <- markers.conserved
        tmp.table <- exp(tmp.table) - 1
        conserved.table$other_mean <- log(rowMeans(as.matrix(tmp.table[markers.conserved, , drop=FALSE])) + 1)

        message("mean levels added")

        print(dim(conserved.table))

        conserved.table <- conserved.table[order(conserved.table$p.adj),]

    } else {
        print("No conserved markers.")
        conserved.table <- data.frame(
            cluster = character(0),
            gene = character(0),
            gene_id = character(0),
            p_val = character(0),
            p.adj = character(0),
            avg_logFC = character(0),
            pct.1 = character(0),
            pct.2 = character(0),
            cluster_mean = character(0),
            other_mean = character(0)
        )
    }

    print("Find(Conserved)Markers complete")
    print(dim(conserved.table))

    if(opt$testfactor=="none"){
        prefix <- "markers.cluster"
    } else {
        prefix <- sprintf("markers.between.%s.cluster", opt$testfactor)
    }

    out_path <- file.path(
            opt$outdir,
            paste(prefix,opt$cluster,"tsv","gz",sep="."))

    print(paste("Saving markers to:", out_path))

    write.table(conserved.table,
                gzfile(out_path),
                quote=F,sep="\t",row.names=F)

    ## write out the ensembl gene_ids of the "universe" for downstream
    ## geneset analysis
    ##
    ## * these are taken to be the genes passing min.pct
    ##   _but not_ min.diff.pct or thresh.use

    message("head of background lists:")
    lapply(background.conserved.list, function(x){print(head(x, 5))})
    # Debug: RDS list of background genes
    # saveRDS(background.conserved.list, file.path(opt$outdir, "background.conserved.list.rds"))

    # intersect of background lists
    background.count.genes <- table(
        unlist(background.conserved.list)) == length(background.conserved.list
    )
    background.conserved.genes <- names(background.count.genes)[background.count.genes]
    message("background.conserved.genes")
    print(length(background.conserved.genes))
    print(tail(background.conserved.genes, 5))

    if(have_gene_ids)
    {
        print("Writing universe gene ids using @misc$gene_id")
        universe_gene_ids = s@misc[background.conserved.genes,"gene_id"]
    } else {
        print("Writing universe gene ids using annotation file")
        ## avoid writing out NA values.
        background.take = background.use[background.conserved.genes %in% rownames(ann)]
        universe_gene_ids = ann[background.take,"ensembl_id"]
    }

    universe <- data.frame(gene_id = universe_gene_ids)

    universe_path <- file.path(
        opt$outdir,
        paste(prefix,opt$cluster,"universe","tsv","gz",sep="."))

    write.table(universe,
                gzfile(universe_path),
                quote=FALSE, sep="\t", row.names=FALSE)
}

message("seurat_FindMarkers.R final default assay: ", DefaultAssay(s))

message("completed")
