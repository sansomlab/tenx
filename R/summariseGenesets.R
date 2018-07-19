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
    make_option(c("--minoddsratio"), type="integer", default=1.5,
                help="The minimum odds ratio."),
    make_option(c("--gmt_names"), default="none",
                help="comma separated list of names for the gmt files"),
    make_option(c("--clustertype"),default="cluster",
                help="will be used e.g. in plot labels"),
    make_option(c("--project"), default="SeuratAnalysis",
                help="project name"),
    make_option(c("--prefix"), default="genesets",
                help="expected prefix for source files"),
    make_option(c("--outfile"), default="none",
                help="outfile")
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
    contents <- NULL

    print(paste("Processing:", geneset,"annotations."))
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
        print(paste("Working on cluster: ", cluster))

        fn = paste0(opt$genesetdir,"/",opt$prefix,".",cluster,".",geneset,".txt.gz")
        if(file.exists(fn))
            {
                temp = read.table(gzfile(fn),sep="\t",header=T,as.is=T,quote="")

                if(nrow(temp)>0)
                {
                    temp$cluster <- cluster
                }
                else {
                    print(paste0("zero rows for cluster: ",cluster))
                }

                if(begin==T)
                {
                    contents <- temp
                    begin <- F
                }
                else {
                    contents <- rbind(contents,temp)
                }
            } else {
                print(paste("Skipping ",fn,"(file not found)",sep="\t"))
            }

    }

    make_plot = FALSE

    if(!is.null(contents))
    {

        ## Filter out genesets we do not wish to consider
        contents <- contents[contents$n_fg >= opt$mingenes,]

        ## Compute adjusted p-values
        ## The input table contains _all tested genesets_
        ## (with n_fg >= opt$mingenes) regardless of p value.
        contents$p.adj <- p.adjust(contents$p.val, method=opt$padjustmethod)

        ## Compute the number of clusters in which each geneset is signficantly enriched

        if(opt$useadjusted)
        {
            nsigtmp <- contents[contents$p.adj < opt$pvaluethreshold,]
        } else {
            nsigtmp <- contents[contents$p.val < opt$pvaluethreshold,]
        }

        id_tab <- table(nsigtmp$geneset_id)

        contents$n_clust_sig <- id_tab[contents$geneset_id]
        contents$n_clust_sig[is.na(contents$n_clust_sig)] <- 0

        ## Sort by p value
        contents <- contents[order(contents$cluster,contents$p.val),]

        ## Tidy up the frame
        firstcols <- c("cluster","geneset_id","description","p.adj","p.val","odds.ratio","n_clust_sig","n_fg","n_bg")
        firstcols <- firstcols[firstcols %in% colnames(contents)]
        othercols <- colnames(contents)[!colnames(contents) %in% firstcols]
        contents <- contents[,c(firstcols,othercols)]

        numeric_cols <- colnames(contents)[sapply(contents, is.numeric)]

        for(numeric_col in numeric_cols)
        {
            ## set to 3 sf
            xx <- contents[[numeric_col]]

            nas <- is.na(xx)
            ints <- all((xx - round(xx)) == 0)
            xx[xx<1000 & !nas & !ints] <- signif(xx[xx<1000 & !nas & !ints],digits=3)
            xx[xx>=1000 & !nas] <- round(xx[xx>=1000 & !nas],digits=0)
            xx[ints] <- as.integer(xx[ints])

            contents[[numeric_col]] <- xx
        }

        ## Add the results to the worksheet
        addWorksheet(wb,geneset)
        setColWidths(wb,geneset,cols=1:ncol(contents),widths=10)
        hs <- createStyle(textDecoration = "BOLD")
        writeData(wb, geneset, contents, withFilter = T, headerStyle=hs)

        ## prepare for writing out a latex summary table (unique pathways)
        ## and geneset heatmaps (can be shared)
        for(clust in unique(as.character(contents$cluster)))
        {

            temp <- contents[contents$cluster==clust,]
            nrows <- nrow(temp)
            if(nrows==0) { next }
            temp <- temp[1:min(nrows,5),]

            if(!"description" %in% colnames(temp))
            {
                temp$description <-temp$geneset_id
            }

            ## trim long descriptions
            maxl <- 45
            xx <- temp$description

            xx <- temp$description
            xx <- formatDescriptions(xx, c("REACTOME_", "BIOCARTA_"), maxl)
            temp$description <- xx

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

        results_table <- contents

        ## catch case where there is nothing to plot.
        if(opt$useadjusted)
        {
            nsig = nrow(results_table[results_table$p.adj < opt$pvaluethreshold &
                                      results_table$odds.ratio >= opt$minoddsratio,])
        } else {
            nsig =  nrow(results_table[results_table$p.val < opt$pvaluethreshold &
                                       results_table$odds.ratio >= opt$minoddsratio,])
        }

        if(nsig > 0) { make_plot <- TRUE }
    }

    plotfn <- gsub("xlsx",geneset,opt$outfile)
    if(make_plot)
    {

        plot_fn <- function()
        {
            sampleEnrichmentHeatmap(results_table,
                                    max_rows=50,
                                    min_genes=opt$mingenes,
                                    p_col="p.val",
                                    adjust_pvalues=opt$useadjusted,
                                    padjust_method=opt$padjustmethod,
                                    pvalue_threshold=opt$pvaluethreshold,
                                    min_odds_ratio=opt$minoddsratio,
                                    maxl=45,
                                    show_common=opt$showcommon,
                                    sample_id_col="cluster",
                                    sample_ids=c(first:last),
                                    title=geneset)
        }

        save_plots(plotfn,
                   plot_fn=plot_fn,
                   width=8,
                   height=8)

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
    tex <- c(tex,getFigureTex(basename(plotfn), caption)) #,plot_dir_var="."))
    tex <- c(tex,"\n")
}

fig_file <- gsub("xlsx","figure.tex",opt$outfile)
writeTex(fig_file,tex)

saveWorkbook(wb, file=opt$outfile, overwrite=T)

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

ltab_file <- gsub("xlsx","table.tex",opt$outfile)

out <- sprintfResults(out)

xtab <- xtable(out, caption="The top (lowest p-value) genesets found (uniquely) in each cluster")

print(xtab,
      include.rownames=F,
      hline.after=hlines,
      file=ltab_file,
      tabular.environment="longtable",
      size="\\fontsize{6pt}{9pt}\\selectfont")
