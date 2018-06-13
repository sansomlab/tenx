## Perform Fisher (hypergeometric) ORA of seurat clusters


# Libraries ----

stopifnot(
    require(optparse),
    require(gsfisher)
)


## deal with the options
option_list <- list(
#    make_option(c("--diffexpall"), default="none",
#                help="table of de analysis of cluster"),
    make_option(c("--markers"), default="none",
                help="the marker summary table"),
    make_option(c("--universe"), default="none",
                help="the list of gene_ids comprising the universe"),
    make_option(c("--species"), default="none",
                help="species: hs or mm"),
    make_option(c("--annotation"),default="none",
                help="entrez_id,ensembl_id,gene_name tab sep"),
    make_option(c("--kegg_pathways"), default="none",
                help="R object containing kegg_pathways list"),
    make_option(c("--gmt_files"), default="none",
                help="comma separated list of gmt files"),
    make_option(c("--gmt_names"), default="none",
                help="comma separated list of names for the gmt files"),
    make_option(c("--cluster"), type="integer", default=0,
                help="the number of the cluster being analysed"),
    make_option(c("--adjpthreshold"),type="double",default=0.1,
                help="p.adj threshold for de genes"),
    make_option(c("--direction"), default="positive",
                  help="direction of changes to consider: positive|negative|both"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--prefix"), default="genesets",
                help="prefix for out files"),
    make_option(c("--outdir"), default="seurat.out.dir",
                help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

print("Running with options:")
print(opt)


## set the run specs
run_specs <- paste(opt$numpcs,opt$resolution,opt$algorithm,opt$testuse,sep="_")
background <- read.table(gzfile(opt$universe), header=T, as.is=T, sep="\t")$gene_id
if(length(background)<100)
{
    stop("Less than 100 genes in the gene universe!!")
}


message("retrieving the gene lists")
de <- read.table(gzfile(opt$markers), header=T, as.is=T, sep="\t")

cluster_de <- de[de$cluster==opt$cluster & de$p.adj<=opt$adjpthreshold,]

## only positive markers considered
if(opt$direction=="positive")
{
    foreground <- unique(cluster_de$gene_id[cluster_de$avg_logFC>0])
} else if (opt$direction=="negative")
{
    foreground <- unique(cluster_de$gene_id[cluster_de$avg_logFC<0])
} else if (opt$direction=="both")
{
    foreground <- unique(cluster_de$gene_id)
} else {
    stop("Direction paramater value not recognised")
}

universe <- unique(c(foreground, background))

print(paste0("no ensembl_ids in foreground: ", length(foreground)))
print(paste0("no ensembl_ids in universe: ", length(universe)))

message("loading annotations")
anno <- read.table(gzfile(opt$annotation), as.is=T, sep="\t", header=T)
print(head(anno))

## get the foreground and universe gene lists
fg_entrez <- anno$entrez_id[anno$ensembl_id %in% foreground]
u_entrez <- anno$entrez_id[anno$ensembl_id %in% universe]

fg_entrez <- as.character(fg_entrez[!is.na(fg_entrez)])
u_entrez <- as.character(u_entrez[!is.na(u_entrez)])

print(paste0("no entrez_ids in foreground: ", length(fg_entrez)))
print(paste0("no entrez_ids in universe: ", length(u_entrez)))

message("performing enrichment tests")
## perform the enrichment tests
if(length(fg_entrez)>0)
    {
        outPrefix = paste0(opt$outdir,"/", opt$prefix, ".",opt$cluster)
        kegg_pathways <- readRDS(opt$kegg_pathways)

        gmts <- list()

        if(!opt$gmt_names == "none" | !opt$gmt_files == "none")
        {
            ## parse the list of gmts

            gmt_names <- strsplit(opt$gmt_names,",")[[1]]
            gmt_files <- strsplit(opt$gmt_files,",")[[1]]

            if(length(gmt_names) != length(gmt_files))
            {
                stop("the same number of gmt_names and gmt_files must be specified")
            }

            for(i in 1:length(gmt_names))
            {
                gmts[[gmt_names[i]]] <- gmt_files[i]
            }
        }

        ## run the analysis.
        results <- analyseGenesets(fg_entrez,
                                   u_entrez,
                                   gene_id_type="entrez",
                                   kegg_pathways=kegg_pathways,
                                   gmt_files=gmts,
                                   species=opt$species)

        ## outPrefix,
        for(geneset in names(results))
            {
                write.table(results[[geneset]],
                            gzfile(paste(outPrefix,geneset,"txt","gz",sep=".")),
                            row.names=FALSE, col.names=TRUE, quote=FALSE,
                            sep="\t")
            }

    } else {
    print("cluster has no significant genes")
    }

message("enrichment tests complete")
