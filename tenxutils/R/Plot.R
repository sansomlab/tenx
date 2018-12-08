## function to return annotated MA plot as ggplot object
#' Make an MA-style plot
#' @param data Matrix or dataframe
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param m_col Column containing M value (log2 ratio)
#' @param a_col Column containing A value (log2 mean expression)
#' @param p_col Column containing p-values.
#' @param top_col Column indicating top rows to highlight.
#' @param label_col Column containing labels for highlighted genes.
plotMA <- function(data, xlab="xlab", ylab="ylab",
                   m_col="M", a_col="A", p_col="p.adj",
                   top_col="top", label_col="gene")
{
    data <- categoriseGenes(data)

    data <- data.frame(data)
    print("====")
    print(head(data))
    top <- data[data[[top_col]] == TRUE & data[[p_col]] < 0.05,]
    print(top)

    midpoint=mean(range(-log10(data[[p_col]][data[[p_col]]>0])))

    gp <- ggplot(data, aes(get(a_col), get(m_col), color=-log10(get(p_col))))
    gp <- gp + geom_point(alpha=1, size=0.75)
    gp <- gp + scale_color_gradient2(low="grey",mid="red",high="black",
                                     midpoint=midpoint, guide=FALSE)
    gp <- gp + geom_text_repel(data=top, aes_string(label=label_col), color="black",
                               min.segment.length=0, size=3)
    gp <- gp + geom_hline(yintercept=c(-1,1), linetype="dashed", color="grey")
    gp <- gp + geom_hline(yintercept=0, linetype="dashed")

    gp <- gp + xlab(xlab) + ylab(ylab)

    gp
}

#' Make an volcano-style plot
#' @param data Matrix or dataframe
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param m_col Column containing M value (log2 ratio)
#' @param a_col Column containing A value (log2 mean expression)
#' @param p_col Column containing p-values.
#' @param top_col Column indicating top rows to highlight.
#' @param label_col Column containing labels for highlighted genes.
plotVolcano <- function(data, xlab="xlab", ylab="ylab",
                        m_col="M", a_col="A", p_col="p.adj",
                        top_col="top", label_col="gene")
{

    data <- categoriseGenes(data)

    data <- data.frame(data)
    top <- data[data[[top_col]] == TRUE & data[[p_col]] < 0.05,]

    midpoint=mean(range(-log10(data[[p_col]][data[[p_col]]>0])))

    gp <- ggplot(data, aes(get(m_col), -log10(get(p_col)), color=-log10(get(p_col))))
    gp <- gp + geom_point(alpha=1, size=0.75)
    gp <- gp + scale_color_gradient2(low="grey",mid="red",high="black",
                                     midpoint=midpoint, guide=FALSE)
    gp <- gp + geom_text_repel(data=top, aes_string(label=label_col), color="black",
                               min.segment.length=0,size=3)
    gp <- gp + geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey")
    gp <- gp + geom_vline(xintercept=c(0), linetype="dashed", color="black")
    gp <- gp + geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey")

    gp <- gp + ylab(ylab) + xlab(xlab)

    gp
}

#' Make a plot of frequency versus expression level
#' @param data Matrix or dataframe
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param m_col Column containing M value (log2 exprs ratio)
#' @param freq_col Column containing L value (log2 freq ratio)
#' @param p_col Column containing p-values.
#' @param label_col Column containing labels for highlighted genes.
plotFvE <- function(data, p_col="p.adj", label_col="gene",
                    m_col="M",freq_col="L",
                    xlab="change in expression (log2)",
                    ylab="change in frequency (log2)")
{
    data <- categoriseGenes(data, use_fc=FALSE, ngenes=12)

    data$scolor <- "1"

    data$scolor[data$sig == TRUE] <- "2"
    data$scolor[data$top == TRUE] <- "3"

    data <- data.frame(data)
    top <- data[data$top == TRUE & data$sig,]

    midpoint=mean(range(-log10(data[[p_col]][data[[p_col]]>0])))

    gp <- ggplot(data, aes(get(freq_col), get(m_col), color=-log10(get(p_col))))
    gp <- gp + geom_point(alpha=1, size=0.75)

    gp <- gp + geom_text_repel(data=top, aes_string(label=label_col),
                               color="black", min.segment.length=0, size=3)

    gp <- gp + geom_hline(yintercept=0, linetype="dashed", color="grey")
    gp <- gp + geom_vline(xintercept=0, linetype="dashed", color="grey")
    gp <- gp + scale_color_gradient2(low="grey",mid="red",high="black",midpoint=midpoint)
    gp <- gp + ylab(xlab)
    gp <- gp + xlab(ylab)

    gp
}


#' Function to make sets of violin plots.
#' @param data The data
#' @param seurat_object The seurat object
#' @param cluster_ids The cluster ids
#' @param type Either "positive" or "negative"
#' @param group.by A factor to group by
#' @param ident.include Identity to include
#' @param ncol Number of columns in the figure
#' @param m_col Column containing M value (log2 exprs ratio)
#' @param use.minfc Use minimum foldchange
#' @param minfc_col Column containing the minimum foldchange
#' @param maxfc_col Column containing the maximum foldchage
#' @param avgfc_col Column containing the average foldchage
#' @param p_col The column containing the p-value
#' @param id_col The column containing unique gene identifiers
plotViolins <- function(data, seurat_object,
                        cluster_ids, type="positive",
                        group.by=NULL, ident.include=NULL,
                        ncol=ncol,
                        use.minfc=FALSE,
                        minfc_col="min_logFC", maxfc_col="max_logFC",
                        avgfc_col="avg_logFC",
                        m_col="avg_logFC", p_col="p.adj",
                        id_col="gene")
{

    ## ensure the gene ids are unique.
    data[[id_col]] <- make.unique(data[[id_col]])

    message("dimensions of data passed to plotViolins()")
    print(dim(data))

    if(type=="positive")
    {
        tmp <- data[data[[m_col]] > 0,]
    } else if (type == "negative")
    {
        tmp <- data[data[[m_col]] < 0,]
    } else {
        stop("type argument to plotViolins() not recognised")
    }

    ## make first sub-figure (4 columns, 3 rows)
    tmp <- tmp[order(tmp[[p_col]]),]

    genes_a <- tmp[[id_col]][min(1,length(tmp[[id_col]])):min(12,length(tmp[[id_col]]))]
    n_a <- length(genes_a)

    message("number of genes for first panel: ", n_a)

    nrow_a <- ceiling(n_a / 4)

    if(n_a > 0)
    {
        message("Calling VlnPlot for top genes by p-value")
        gpa <- VlnPlot(object=seurat_object, features.plot=genes_a,
                       size.x.use=6,size.y.use=6,size.title.use=10, point.size.use=0.1,
                       nCol=ncol,
                       y.log=TRUE, do.return=TRUE,
                       group.by=group.by, ident.include=ident.include)

        gpa_exists <- TRUE

    } else {
        gpa_exists <- FALSE
        gpa <- FALSE
        nrow_a <- 0
    }

    ## if n_a is less than 12, there are no more genes to show.
    if(n_a == 12)
    {

        ## make second sub-figure
        tmp <- tmp[!tmp[[id_col]] %in% genes_a,]

        ## order by largest fold change
        ## tmp <- tmp[rev(order(abs(tmp$avg_logFC))),]

        ## subset to positive or negative markers
        if(type == "positive")
        {
            ## order the positive markers by fold change
            if(use.minfc)
            {
                tmp <- tmp[tmp[[minfc_col]] > 0,]
                tmp <- tmp[rev(order(tmp[[minfc_col]])),]
            } else {
                tmp <- tmp[rev(order(tmp[[avgfc_col]])),]
            }
        } else {

            ## order the negative markers by fold change
            if(use.minfc)
            {
                tmp <- tmp[tmp[[maxfc_col]] < 0,]
                tmp <- tmp[order(tmp[[maxfc_col]]),]
            } else {
                tmp <- tmp[order(tmp[[avgfc_col]]),]
            }
        }

        genes_b <- tmp[[id_col]][min(1,length(tmp[[id_col]])):min(12,length(tmp[[id_col]]))]
        n_b <- length(genes_b)
        nrow_b <- ceiling(n_b/4)

        message("number of genes for second panel: ", n_b)

        if(n_b > 0)
        {

            message("Calling VlnPlot for top genes by fold change")
            gpb <- VlnPlot(object=seurat_object, features.plot=genes_b,
                           size.x.use=6,size.y.use=6,size.title.use=10, point.size.use=0.1,
                           nCol=ncol,
                           y.log=TRUE, do.return=TRUE,
                           group.by=group.by, ident.include= ident.include)

            gpb_exists <- TRUE

        } else {
            gpb <- FALSE
            gpb_exists <- FALSE
            nrow_b <- 0
        }
    } else {
        gpb <- FALSE
        gpb_exists <- FALSE
        nrow_b <- 0

    }

    return(list(gpa=gpa, nrow_a=nrow_a, gpa_exists=gpa_exists,
                gpb=gpb, nrow_b=nrow_b, gpb_exists=gpb_exists))
}



#' Make a latex section containing a set of violin plots
#' @param data The data
#' @param seurat_object The seurat object
#' @param cluster_ids The cluster ids
#' @param type Either "positive" or "negative"
#' @param group.by A factor to group by
#' @param ident.include Identity to include
#' @param ncol Number of columns in the figure
#' @param use.minfc Use minimum foldchange
#' @param outdir The directory for the output
#' @param analysis_title Title for this section
#' @param fc_type The metric by which the violin plots are ordered
violinPlotSection <- function(data, seurat_object, cluster_ids, type="positive",
                              group.by=opt$testfactor,
                              ident.include=opt$identinclude, ncol=ncol,
                              outdir=opt$outdir,
                              analysis_title="violin plots", fc_type="fold change",
                              plot_dir_var="plotsDir",
                              use.minfc=FALSE)
{

    ## get the violin plots
    violin_plots <- plotViolins(data, seurat_object, cluster_ids, type=type,
                                group.by=group.by,
                                ident.include=ident.include, ncol=ncol,
                                use.minfc=use.minfc)


    tex <- c()

    subsectionTitle <- getSubsectionTex(paste0("Cluster ",
                                               cluster, " violin plots: ",
                                               type,
                                               analysis_title))
    tex <- c(tex, subsectionTitle)

    ## make the figure
    tex <- c(tex, "\\begin{figure}[H]")

    if(violin_plots$gpa_exists)
    {

        ## first subfigure
        violin_fn <- gsub("type", paste0(type,".padj"), violin_fn_template)
        violin_path <- file.path(outdir, violin_fn )

        h <- max(violin_plots$nrow_a * 2, 3)

        save_ggplots(violin_path,
                     violin_plots$gpa,
                     width=10,
                     height=h)

        caption <- paste0("Top ", type, analysis_title,
                          " ordered by p-value, cluster: ",cluster)

        tex <- c(tex, getSubFigureTex(violin_fn, caption, plot_dir_var=plot_dir_var))

        if(violin_plots$gpb_exists)
        {

            ## second subfigure
            violin_fn <- gsub("type",paste0(type,".fc"), violin_fn_template)
            violin_path <- file.path(opt$outdir, violin_fn)

            h <- max(violin_plots$nrow_b * 2, 3)

            save_ggplots(violin_path,
                         violin_plots$gpb,
                         width=10,
                         height=h)

            caption <- paste0("Additional ", type, analysis_title, " ordered by ", fc_type, ", cluster: ",cluster)
            tex <- c(tex, getSubFigureTex(violin_fn, caption, plot_dir_var=plot_dir_var))
        }

    } else {
        tex <- c(tex, "No significant genes")
    }

    tex <- c(tex, "\\end{figure}")
    tex

}

#' Extract a legend from a ggplot
#' @param a.ggplot A ggplot obect
g_legend<-function(a.ggplot){
    tmp <- ggplot_gtable(ggplot_build(a.ggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

#' function to return colors matching ggplot2
#' see: https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @param n Number of colours to return.
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


## save png and pdf version of base graphics plots
#' Wrapper function for saving base R plots in both pdf and png format
#' @param filepath The name of the file to write to (without .pdf or .png extension)
#' @param plot_fn The plotting function to call.
#' @param width The width of the plots.
#' @param height The height of the plots.
#' @param res The plot resolution.
#' @param to_pdf Save pdf version of plot (TRUE|FALSE).
#' @param to_png Save png version of plot (TRUE|FALSE).
save_plots <- function(filepath=NULL,
                       plot_fn=NULL,
                       width=6,
                       height=8,
                       res=300,
                       to_pdf=TRUE,
                       to_png=TRUE)
{

    ## save the pdf plot
    if(to_pdf)
    {
        cairo_pdf(paste(filepath,"pdf",sep="."),
                  width=width,
                  height=height)
        do.call(plot_fn,list())
        dev.off()

    }

    ## save the png plot
    if(to_png)
    {
        png(paste(filepath,"png",sep="."),
            width=width,
            height=height,
            units="in",
            res=res)
        do.call(plot_fn,list())
        dev.off()
    }

}


## save png and pdf versions of ggplot plots
## save png and pdf version of base graphics plots
#' Wrapper function for saving base R plots in both pdf and png format
#' @param filepath The name of the file to write to (without .pdf or .png extension)
#' @param gp The ggplot object.
#' @param width The width of the plots.
#' @param height The height of the plots.
#' @param dpi The plot resolution.
#' @param to_pdf Save pdf version of plot (TRUE|FALSE).
#' @param to_png Save png version of plot (TRUE|FALSE).
save_ggplots <- function(filepath=NULL,
                         gp=NULL,
                         width=6,
                         height=8,
                         dpi=300,
                         to_pdf=TRUE,
                         to_png=TRUE)
{
    if(to_png)
    {
        ggsave(paste(filepath,"png",sep="."),
               gp,
               device='png',
               width=width,
               height=height,
               units="in",
               dpi=dpi)
    }

    if(to_pdf)
    {
        ggsave(paste(filepath,"pdf",sep="."),
               gp,
               device=cairo_pdf,
               width=width,
               height=height,
               units="in",
               dpi=dpi)
    }
}

#' Plot
#'
#' @param matrixUMI
#' @param metadata
#' @param basename
#'
#' @return
#' @export
#'
#' @examples
plotDownsampling <- function(matrixUMI, metadata, basename) {

    # Total UMI count per cell
    nUMIs <- Matrix::colSums(matrixUMI)

    # Collate total UMI and cell rank for each cell of each sample
    inputStats <- do.call("rbind", lapply(
        levels(metadata$sample_id),
        function(id){
            codes_id <- subset(metadata, sample_id == id, "barcode", drop=TRUE)
            nUMIs_id <- nUMIs[codes_id]
            data.frame(
                nUMIs=sort(nUMIs_id, decreasing=TRUE),
                CellRank=seq_along(nUMIs_id),
                sample=id
            )
        }))

    # UMI vs Rank ----

    gg <- ggplot(inputStats) +
        geom_line(aes(CellRank, nUMIs, colour=sample), size=0.25) +
        scale_x_log10(
            limits=c(1, max(inputStats$CellRank)) #, breaks=c(1,10,1E2,1E3,1E4,1E5,1E6)
        ) +
        scale_y_log10(
            limits=c(1, 10^ceiling(log10(max(nUMIs)))) #, breaks=c(1,2,5,10,20,50,1E2,2E2,5E2,1E3,2E3,5E3,1E4,2E4,5E4,1E5,2E5,5E5,1E6)
        ) +
        annotation_logticks() +
        labs(y="UMI count", x="Cell rank") +
        theme_bw() +
        theme(
            # axis.text=element_text(size=rel(0.5)),
            # legend.title=element_text(size=rel(0.75)),
            # legend.text=element_text(size=rel(0.75)),
            panel.grid.major=element_line(size=0.1, color="grey"),
            panel.grid.minor=element_blank()
        )

    ggsave(
        file.path(opt$outdir, sprintf("%s_ranked.pdf", basename)), gg,
        width=7, height=5
    )

    # UMI violion plot ----

    gg <- ggplot(inputStats) +
        geom_violin(
            aes(sample, nUMIs, colour=sample),
            draw_quantiles=c(0.5), size=0.25) +
        scale_y_log10(
            limits=c(1, 10^ceiling(log10(max(nUMIs)))) #, breaks=c(1,2,5,10,20,50,1E2,2E2,5E2,1E3,2E3,5E3,1E4,2E4,5E4,1E5,2E5,5E5,1E6)
        ) +
        annotation_logticks(sides="l") +
        labs(y="UMI count", x="Sample") +
        theme_bw() +
        theme(
            # axis.text=element_text(size=rel(0.5)),
            # legend.title=element_text(size=rel(0.75)),
            # legend.text=element_text(size=rel(0.75)),
            panel.grid.major=element_line(size=0.1, color="grey"),
            panel.grid.minor=element_blank(),
            axis.text.x=element_text(angle=90)
        )

    ggsave(
        file.path(opt$outdir, sprintf("%s_violin.pdf", basename)), gg,
        width=7, height=5
    )

    return(TRUE)
}
