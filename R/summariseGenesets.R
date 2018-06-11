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
    make_option(c("--pthreshold"),type="double",default=0.05,
                help="raw p threshold filtering sets"),
    make_option(c("--mingenes"), type="integer", default=2,
                help="min no. genes in foreground set"),
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

## TODO: Detect automatically
genesets = c("GO.BP","GO.MF","GO.CC",
             "KEGG",
              strsplit(opt$gmt_names,",")[[1]])

## set up workbook.
wb <- createWorkbook()

ltabs <- list()
hmaps <- list()
tex <- c()

for(geneset in genesets)
{
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

    if(nrow(contents) == 0)
    {
        print(paste0("Warning: no entries for geneset: ",geneset))
        next
    }

    ## Filter the genesets
    use_adjust_pvalues = FALSE
    pvalue_threshold = 0.1
    show_common = TRUE
    padjust_method="BH"
    total_n_clust <- length(unique(contents$cluster))

    contents <- contents[contents$n_fg >= opt$mingenes,]

    contents$p.adj <- p.adjust(contents$p.val, method=padjust_method)

    if(use_adjust_pvalues)
    {
      ## compute FDR accross all samples
      contents <- contents[contents$p.adj < pvalue_threshold
                                     & !is.na(contents$p.adj),]
    } else {
      contents <- contents[contents$p.val < pvalue_threshold
                                     & !is.na(contents$p.val),]
    }

    if(!show_common)
    {
      contents <- contents[contents$n_sample < total_n_sample,]
    }


    ## Compute the number of clusters in which the geneset is enriched
    id_tab <- table(contents$geneset_id)
    contents$n_clust <- id_tab[contents$geneset_id]

    ## Sort by p value
    contents <- contents[order(contents$cluster,contents$p.val),]

    ## Tidy up the frame
    firstcols <- c("cluster","geneset_id","description","p.adj","p.val","odds.ratio","n_clust","n_fg","n_bg")
    firstcols <- firstcols[firstcols %in% colnames(contents)]
    othercols <- colnames(contents)[!colnames(contents) %in% firstcols]
    contents <- contents[,c(firstcols,othercols)]

    numeric_cols <- colnames(contents)[sapply(contents, is.numeric)]

    for(numeric_col in numeric_cols)
    {
        ## set to 3 sf
        xx <- contents[[numeric_col]]

        nas <- is.na(xx)
        xx[xx<1000 & !nas] <- signif(xx[xx<1000 & !nas],digits=3)
        xx[xx>=1000 & !nas] <- round(xx[xx>=1000 & !nas],digits=0)

        contents[[numeric_col]] <- xx
    }

    addWorksheet(wb,geneset)
    setColWidths(wb,geneset,cols=1:ncol(contents),widths=10)
    hs <- createStyle(textDecoration = "BOLD")
    writeData(wb, geneset, contents, withFilter = T, headerStyle=hs)

    ## prepare for writing out a latex summary table (unique pathways)
    ## and geneset heatmaps (can be shared)
    for(clust in unique(as.character(contents$cluster)))
    {

        # deal with latex
        temp <- contents[contents$cluster==clust & contents$n_clust==1 & contents$n_fg >=3,]
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
        xx[is.na(xx)] <- "n/a"
        xx[nchar(xx)>maxl] <- paste0(strtrim(xx[nchar(xx)>maxl],maxl),"...")
        temp$description <- xx


        temp_names <- colnames(temp)
        temp$type <- geneset
        temp <- temp[,c("type",temp_names)]

        temp <- temp[,c("type","description","p.val","p.adj","odds.ratio","n_fg")]

        if(clust %in% names(ltabs))
        {
            ltabs[[clust]] <- rbind(ltabs[[clust]],temp)
        }
        else {

            ltabs[[clust]] <- temp
        }
    }

    results_table <- contents

    plotfn <- gsub("xlsx",geneset,opt$outfile)

    print(dim(results_table))
    print(head(results_table))

    if(nrow(results_table) > 0)
    {

    plot_fn <- function()
    {
      sampleEnrichmentHeatmap(results_table,
                              max_rows=50,
                              min_genes=opt$mingenes,
                              adjust_pvalues=use_adjust_pvalues,
                              padjust_method=padjust_method,
                              pvalue_threshold=pvalue_threshold,
                              maxl=45,
                              show_common=show_common,
                              sample_id_col="cluster",
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
            text(0.5,0.5,paste0("no data for:\n",geneset))
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
xtab <- xtable(out, caption="The top (lowest p-value) genesets found (uniquely) in each cluster")
print(xtab,
      include.rownames=F,
      hline.after=hlines,
      file=ltab_file,
      tabular.environment="longtable",
      size="\\fontsize{6pt}{9pt}\\selectfont")
