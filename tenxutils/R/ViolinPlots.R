#' Get data for a set of violin plots
#' @param seurat_object A seurat object
#' @param metadata A vector of metadata columns to retrieve
#' 
#' @export
#' 
getViolinData <- function(seurat_object, 
                          genes, 
                          metadata=NULL,
                          clusters=NULL)
{
  require(reshape2)

  genes <-genes[genes %in% rownames(seurat_object@data)]
    
  data <- as.data.frame(as.matrix(seurat_object@data[genes,]))

  data$gene <- as.vector(rownames(data))
  
  ggData <- melt(data,id.vars="gene")
  ggData$cluster <- as.numeric(as.vector(seurat_object@ident[ggData$variable]))
  
  if(!is.null(clusters))
  {
    ggData <- ggData[ggData$cluster %in% clusters,]
  }
  
  if(length(metadata)>0)
  {
    ggData[,metadata] <- seurat_object@meta.data[ggData$variable, metadata]
  }
  ggData
}


#' Remove spacer grobs from a facetted ggplot plot.
#' @param gp The ggplot object containing the violin plots
#' @param ncol The number of columns (as was passed to facet wrap)
#' @param nmissing  The number of facet panels missing
#' 
#' @export
#' 
removeGrobs <- function(gpGrob, ncol, nmissing)
{
  
  to_remove <- c()
  for(i in (1+(ncol-nmissing)):ncol)
  {
    to_remove <- c(to_remove, paste("panel",i,"1",sep="-"),
                   paste("strip-t",i,"1",sep="-"))
  }
  
  # get the grobs that must be removed
  rm_grobs <- gpGrob$layout$name %in% to_remove
  
  # remove grobs
  gpGrob$grobs[rm_grobs] <- NULL
  gpGrob$layout <- gpGrob$layout[!rm_grobs, ]
  
  gpGrob
}

#' Function to draw a block of horizontal violin lots
#' 
#' The number of columns will be respected even if
#' there are fewer genes than columns
#' 
#' @param ggData melted dataframe ready for plotting
#' @param title A title for this set of genes
#' @param ncol  Number of columns, will be passed to facet_wrap
#' @param group A column on which to group plots by
#' @param colors A vector of (fill) colors, one per cluster.
#' @param xlab The label for the x axis.
#' 
#' @export
#' 
makeViolins <- function(ggData, title=NULL, ncol=8, group=NULL,
                        colors=NULL,
                        xlab=NULL)
{
  require(ggplot2)
  require(ggstance)
  theme_set(theme_classic(base_size = 8))
  
  ggData$gene <- factor(ggData$gene)
  
  cluster_levels <- unique(ggData$cluster)
  ggData$cluster <- factor(ggData$cluster,
                           levels=cluster_levels[rev(order(as.numeric(cluster_levels)))])
  
  
  nl  <- length(levels(ggData$gene))
  nrow <- ceiling(nl/ncol)
  
  if(nl < ncol) { to_add <- (ncol - nl) } else { to_add <- 0 }
  
  if(to_add > 0)
  {
    ggData$gene <- factor(ggData$gene, 
                          levels = c(levels(ggData$gene),
                                     paste0(".",c(1:to_add))))
  }
  
  ticks <- function() {function(limits) c(round(min(limits)),max(1,floor(max(limits))))}
  
  if(is.null(group))
  {
    gp <- ggplot(ggData, aes(value, cluster, fill=cluster)) 
  } else {
    ggData[[group]] <- factor(ggData[[group]])
    gp <- ggplot(ggData, aes_string("value", "cluster", fill=group)) 
  }
  gp <- gp + geom_violinh(scale = "width", trim = TRUE) 
  gp <- gp + facet_wrap(~gene, scales="free_x", ncol=ncol, drop=F)
  if(!is.null(title)) {
    gp <- gp + ggtitle(title)
  }
  if(!is.null(colors)) {
    gp <- gp + scale_fill_manual(values = colors)
  }
  gp <- gp + theme_classic()
  gp <- gp + scale_x_continuous(breaks = ticks())
  gp <- gp + xlab(xlab)
  gp <- gp + theme(legend.position = "none", 
                   strip.background = element_rect(fill="grey94",color="grey94"),
                   panel.spacing.x=unit(2, "mm"),
                   axis.line.y = element_blank())
  #  axis.title.x = element_blank())
  gp <- gp + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)
  
  gpGrob <- ggplotGrob(gp) 
  
  if(nl < ncol) {
    gpGrob <- removeGrobs(gpGrob, ncol, ncol - nl)
  } 
  
  gpGrob
}

#' Function to draw grob on a new page
#' 
#' @param grob A grob.
#' 
#' @export
#' 
plotGrob <- function(ggGrob)
{
  require(grid)
  grid.newpage()
  grid.draw(ggGrob)
}


#' @param seurat_object A seurat object
#' @param ggData melted dataframe ready for plotting
#' @param title A title for this set of genes
#' @param ncol  Number of columns, will be passed to facet_wrap
#' @param group A column on which to group plots by
#' @param xlab The label for the x axis.
#' @param colors A vector of (fill) colors, one per cluster.
#' 
#' @export
#' 
plotHorizontalViolins <- function(seurat_object, genes, clusters=NULL, 
                        title=NULL, ncol=8,  
                        xlab="normalised expression level",
                        group=NULL,  colors=NULL, plot=TRUE)

{
  ggData <- getViolinData(seurat_object, 
                          genes=genes, 
                          clusters=clusters,
                          metadata=group)
  
  ggGrob <- makeViolins(ggData, 
                        title=title, 
                        ncol=ncol, 
                        xlab=xlab,
                        group=group,
                        colors=colors)

  if(plot)
  {
    plotGrob(ggGrob)
  } else {
    ggGrob
  }
}
