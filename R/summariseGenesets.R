# Libraries ----

stopifnot(
  require(optparse),
  require(reshape2),
  require(gplots),
  require(openxlsx),
  require(xtable),
  require(tenxutils),
  require(gsfisher)
)

# Options ----

option_list <- list(
    make_option(c("--genesetdir"), default="none",
                help="directory containing the genesets to aggregate"),
    make_option(c("--nclusters"), type="integer", default=0,
                help="the number of the clusters being analysed"),
    make_option(c("--firstcluster"), type="integer", default=0,
                help="clusters might not be zero based..."),
    make_option(c("--pvaluethreshold"),type="double",default=0.05,
                help="p value threshold for filtering sets"),
    make_option(c("--padjustmethod"), default="BH",
                help="The given method is passed to p.adjust"),
    make_option(c("--useadjusted"), default=TRUE,
                help="should adjusted p-values be used for the summary heatmap"),
    make_option(c("--showcommon"), default=TRUE,
                help=paste("Should genesets significantly enriched in all clusters",
                           "be shown in the summary heatmap")),
    make_option(c("--mingenes"), type="integer", default=2,
                help="min no. genes in foreground set"),
    make_option(c("--maxgenes"), type="integer", default=500,
                help="the maximum number of genes allowed per geneset"),
    make_option(c("--minoddsratio"), type="double", default=1.5,
                help="The minimum odds ratio."),
    make_option(c("--gmt_names"), default="none",
                help="comma separated list of names for the gmt files"),
    make_option(c("--show_detailed"), default="none",
                help=paste("comma separated list of names for which to make individual",
                "per-sample/cluster plots")),
    make_option(c("--clustertype"),default="cluster",
                help="will be used e.g. in plot labels"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--prefix"), default="genesets",
                help="expected prefix for source files"),
    make_option(c("--plotdirvar"), default="clusterGenesetsDir",
                help="latex var containing name of the directory with the plots"),
    make_option(c("--outprefix"), default="none",
                help="prefix for outfiles")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


## aggregate geneset types by worksheet
## all clusters in each worksheet

if(opt$gmt_names != "none")
{
    gmt_names <- strsplit(opt$gmt_names,",")[[1]]
} else {
    gmt_names <- c()
}

if(opt$show_detailed != "none")
{
    show_detailed <- strsplit(opt$show_detailed,",")[[1]]
} else {
    show_detailed <- c()
}



## TODO: Detect automatically
genesets = c("GO.BP","GO.MF","GO.CC",
             "KEGG",
              gmt_names)

## set up workbook.
wb <- createWorkbook()

ltabs <- list()
hmaps <- list()
tex <- c()

for(geneset in genesets)
{
    genesets <- NULL

    message(paste("Processing:", geneset,"annotations."))
    begin=T

    if(opt$firstcluster==0)
    {
        first <- 0
        last <- opt$nclusters - 1
        }
    else {
            first <- opt$firstcluster
            last <- opt$nclusters}

    ## build a single table containing the results of the geneset tests for
    ## all clusters
    for(cluster in first:last)
    {
        message(paste("Working on cluster: ", cluster))

        fn = paste0(opt$genesetdir,"/",opt$prefix,".",cluster,".",geneset,".txt.gz")
        if(file.exists(fn))
            {
                temp = read.table(gzfile(fn),sep="\t",header=T,as.is=T,quote="")

                if(nrow(temp)>0)
                {
                    temp$cluster <- cluster
                }
                else {
                    message(paste0("zero rows for cluster: ",cluster))
                }

                if(begin==T)
                {
                    genesets <- temp
                    begin <- F
                }
                else {
                    genesets <- rbind(genesets,temp)
                }
            } else {
                message(paste("Skipping ",fn,"(file not found)",sep="\t"))
            }

    }

    make_plot = FALSE

    if(!is.null(genesets))
        {
            ## Filter out genesets we do not wish to consider
            filtered_genesets <- filterGenesets(genesets,
                                   min_foreground_genes = opt$mingenes,
                                   max_genes_geneset = opt$maxgenes,
                                   min_odds_ratio = opt$minoddsratio,
                                   padjust_method=opt$padjustmethod,
                                   use_adjusted_pvalues=opt$useadjusted,
                                   pvalue_threshold=opt$pvaluethreshold)


            results_table <- filtered_genesets
            } else { results_table <- NULL }

    if(!is.null(results_table) && nrow(results_table) > 0)
    {

        id_tab <- table(results_table$geneset_id)

        results_table$n_clust_sig <- id_tab[results_table$geneset_id]
        results_table$n_clust_sig[is.na(results_table$n_clust_sig)] <- 0

        ## Sort by p value
        results_table <- results_table[order(results_table$cluster,results_table$p.val),]

        ## Tidy up the frame
        firstcols <- c("cluster","geneset_id","description",
                       "p.adj","p.val",
                       "odds.ratio",
                       "n_clust_sig","n_fg","n_bg")

        firstcols <- firstcols[firstcols %in% colnames(results_table)]
        othercols <- colnames(results_table)[!colnames(results_table) %in% firstcols]
        results_table <- results_table[,c(firstcols,othercols)]

        numeric_cols <- colnames(results_table)[sapply(results_table, is.numeric)]

        for(numeric_col in numeric_cols)
        {
            ## set to 3 sf
            xx <- results_table[[numeric_col]]

            nas <- is.na(xx)
            if(any(abs(xx)==Inf))
            {
                ints <- FALSE
            } else {
                ints <- all((xx - round(xx)) == 0)
            }
            xx[xx<1000 & !nas & !ints] <- signif(xx[xx<1000 & !nas & !ints],digits=3)
            xx[xx>=1000 & !nas] <- round(xx[xx>=1000 & !nas],digits=0)

            xx[ints] <- as.integer(xx[ints])

            results_table[[numeric_col]] <- xx
        }

        ## Add the results to the worksheet
        addWorksheet(wb,geneset)
        setColWidths(wb,geneset,cols=1:ncol(results_table),widths=10)
        hs <- createStyle(textDecoration = "BOLD")
        writeData(wb, geneset, results_table, withFilter = T, headerStyle=hs)

        ## prepare for writing out a latex summary table (unique pathways)
        ## and geneset heatmaps (can be shared)
        for(clust in unique(as.character(results_table$cluster)))
        {

            temp <- results_table[results_table$cluster==clust,]
            nrows <- nrow(temp)
            if(nrows==0) { next }
            temp <- temp[1:min(nrows,5),]

            if(!"description" %in% colnames(temp))
            {
                temp$description <-temp$geneset_id
            }

            ## trim long descriptions
            maxl <- 45

            temp$description <- formatDescriptions(temp$description,
                                                   c("REACTOME_", "BIOCARTA_"),
                                                   maxl)

            temp_names <- colnames(temp)
            temp$type <- geneset
            temp <- temp[,c("type",temp_names)]

            temp <- temp[,c("type","description","p.val","p.adj",
                            "n_fg","odds.ratio","n_clust_sig")]

            colnames(temp) <- c("type","description","p.val","p.adj",
                                "n_fg","odds.ratio","n.clust")

            if(clust %in% names(ltabs))
            {
                ltabs[[clust]] <- rbind(ltabs[[clust]],temp)
            }
            else {
                ltabs[[clust]] <- temp
            }
        }

        if(nrow(results_table) > 0) { make_plot <- TRUE }
    }

    plotfn <- paste(opt$outprefix, geneset, sep=".")
    if(make_plot)
    {

        xx <- filtered_genesets

        if(!opt$showcommon)
        {
           tmp <- table(xx$geneset_id)
           xx <- xx[!xx$geneset_id %in% names(tmp)[tmp==opt$nclusters],]
        }

        xx$score <- -log10(xx$p.adj) * log2(xx$odds.ratio)

        genesets_to_show <- getSampleGenesets(xx,
                                              sort_by = "score",
                                              max_rows = 50)

        # add back adjusted p values
        genesets$p.adj <- 1
        genesets[rownames(filtered_genesets),"p.adj"] <- filtered_genesets$p.adj

        message("making sample enrichment dotplot with n=",nrow(genesets)," genesets")

        gp <- sampleEnrichmentDotplot(genesets,
                                      selected_genesets = genesets_to_show,
                                      selection_col = "geneset_id",
                                      sample_levels =c(first:last),
                                      min_dot_size =1, max_dot_size = 6,
                                      maxl = 45,
                                      pvalue_threshold = opt$pvaluethreshold,
                                      title=geneset)

        print(plotfn)
        save_ggplots(plotfn,
                     gp,
                     width=8,
                     height=8)


        message("saved sample enrichement dotplot")

        per_sample_tex = c()
        if(geneset %in% show_detailed)
        {

                ## make the per sample plots
                for(cluster in unique(xx$cluster))
                {
                    tmp <- xx[xx$cluster==cluster,]
                    tmp <- tmp[rev(order(tmp$score)),]

                    max_n_cat = 150

                    if(nrow(tmp)> max_n_cat) { tmp <- tmp[1:max_n_cat,] }

                    if("description" %in% colnames(tmp))
                    {
                       name_col <- "description"
                    } else { name_col <- "geneset_id" }

                    if(nrow(tmp) > 1)
                        {
                        gp <- visualiseClusteredGenesets(tmp,
                                                         highlight=genesets_to_show[genesets_to_show %in% tmp$geneset_id],
                                                         name_col=name_col)

                        detailed_plotfn <- paste(opt$outprefix,
                                                 geneset, "circle_plot", cluster, sep=".")

                        save_ggplots(detailed_plotfn,
                                     gp,
                                     width=10, height=10)

                        caption <- paste("Cluster", cluster, geneset,
                                         "genesets clustered by similarity between over-represented genes.", sep=" ")

                        per_sample_tex <- c(per_sample_tex,
                                            getFigureTex(basename(detailed_plotfn),
                                                         caption,
                                                         plot_dir_var=opt$plotdirvar))

                    }}
        }

    } else {
            # draw an empty plot with an error message
            pngfn <- paste(plotfn, "png", sep=".")
            png(pngfn,width=8,height=8,units="in",res=100)
            plot.new()
            text(0.5,0.5,paste0("no significant genesets for:\n",geneset))
            dev.off()
    }

    caption <- paste("Heatmap of the top", geneset, "genesets", sep=" ")
    tex <- c(tex,getSubsectionTex(geneset))

    tex <- c(tex,getFigureTex(basename(plotfn), caption,
                              plot_dir_var=opt$plotdirvar))

    tex <- c(tex, "\n",
             per_sample_tex, "\n")

}

fig_file <- paste(opt$outprefix,"figure.tex", sep=".")
writeTex(fig_file,tex)

saveWorkbook(wb,
             file=paste(opt$outprefix, "xlsx", sep="."),
             overwrite=T)

begin=T
hlines <- c()


for(cluster in names(ltabs))
{
    temp <- ltabs[[cluster]]
    temp_names <- colnames(temp)
    temp$cluster <- cluster
    temp <- temp[,c("cluster",temp_names)]
    if(begin==T)
    {
        out <- temp
        r <- nrow(temp)
        hlines <- r
        begin <- F
    }
    else {
        out <- rbind(out, temp)
        r <- r + nrow(temp)
        hlines <- c(hlines,r)
         }


}

ltab_file <- paste(opt$outprefix,"table.tex", sep=".")

if(!exists("out"))
{
    out <- data.frame(x=c("no significantly enriched genesets found"))
} else {
    out <- sprintfResults(out)
}

xtab <- xtable(out, caption="The top (lowest p-value) genesets found (uniquely) in each cluster")

print(xtab,
      include.rownames=F,
      hline.after=hlines,
      file=ltab_file,
      tabular.environment="longtable",
      size="\\fontsize{6pt}{9pt}\\selectfont")
