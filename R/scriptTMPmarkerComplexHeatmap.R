message("Loading temp modified function that allows for number of plotted genes to be less than the no of markers.
Function temp is used only if the default (tenexutils::markerComplexHeatmap) fails")
TMPmarkerComplexHeatmap  <-function(seurat_object,
                                    marker_table=NULL,
                                    n_markers=20,
                                    cells_use=NULL,
                                    row_names_gp=10,
                                    sub_group=NULL,
                                    disp_min=-2.5,
                                    disp_max=2.5)
{
  
  # we need the recently added vertical split functionality of ComplexHeatmap.
  if(!packageVersion("ComplexHeatmap")>"1.99")
  {
    stop(paste("ComplexHeatmap version 1.99 or later is needed.",
               "it is avaliable via devtools here: https://github.com/jokergoo/ComplexHeatmap"))
  }
  require(circlize)
  require(RColorBrewer)
  require(colormap)
  top_markers <- marker_table %>% group_by(cluster) %>% top_n(n=n_markers,wt=avg_logFC)
  
  if(is.null(cells_use)) {cell_use <- colnames(GetAssayData(s, slot="scale.data")) }
  
  genes_use <- top_markers$gene[top_markers$gene %in% rownames(GetAssayData(seurat_object,
                                                                            slot="scale.data"))]
  # this check makes sure that the genes fetched from the s object are those scaled, but doesn't
  # align the top_markers dimensions to the same criteria.
  # now, if the full gexp matrix was scaled this should not be a problem, but if you 
  # have scaled only the top var genes for example. that's where the no of rows of top_markers
  # have to match with genes_use.
  # the best scenario, however, is that you scale the entire gexp matrix.
  if(any(!top_markers$gene %in% genes_use) ){ 
    message("subsetting top_markers to the number of genes that were originally scaled, 
            consider re-scaling your object  to include all genes in the scale.data slot")
    top_markers <- top_markers %>% filter(gene %in% genes_use)}
  
  cells_use <- cell_use %in% colnames(GetAssayData(seurat_object, slot="scale.data"))
  
  x <- as.matrix(GetAssayData(seurat_object, slot="scale.data")[genes_use, cells_use])
  
  x <- MinMax(x, min = disp_min, max = disp_max)
  
  clusters <- Idents(seurat_object)[cells_use]
  
  
  # set up the cluster color palette
  nclust <- length(unique(clusters))
  cluster_cols <- gg_color_hue(nclust)
  names(cluster_cols) <- sort(as.numeric(as.character(unique(clusters))))
  
  # get the vector of per-cell cluster names
  cell_clusters <- clusters[colnames(x)]
  
  clusterAnnotation = rowAnnotation(df = data.frame(cluster=top_markers$cluster),
                                    col = list(cluster=cluster_cols),
                                    show_annotation_name = FALSE,
                                    show_legend=FALSE,
                                    width=unit(2,"mm"))
  
  
  
  
  if(!is.null(sub_group))
  {
    
    if(!sub_group %in% colnames(s[[]]))
    {
      stop("specified sub group not found in the metadata")
    }
    
    # set up the subgroup colour palette
    cell_sub_groups <- s[[]][cells_use, sub_group]
    sub_groups <- unique(cell_sub_groups)
    if(length(sub_groups>6)){
    sub_group_cols <- colormap(colormap = colormaps$portland,nshades = length(sub_groups))
    }else{
    sub_group_cols <- brewer.pal(length(sub_groups),"Greys")
    }
    # because the Grey palette returns a minimum of 3 colors..
    sub_group_cols <- sub_group_cols[1:length(sub_groups)]
    names(sub_group_cols) <- sub_groups
    
    print(sub_group_cols)
    # get an ordering vecotor
    cell_order <- order(cell_clusters, cell_sub_groups)
    
    subgroupAnnotation = HeatmapAnnotation(df=data.frame(subgroup=cell_sub_groups[cell_order]),
                                           col = list(subgroup=sub_group_cols),
                                           show_annotation_name = FALSE)
  } else {
    cell_order <- order(cell_clusters)
    subgroupAnnotation <- NULL
  }
  
  x <- x[,cell_order]
  
  # match the Seurat colors
  exprs_cols <- colorRamp2(c(-2.5,0,2.5), c("#FF00FF", "#000000", "#FFFF00"))
  
  ## fudge the row font size to something sensible.
  n_rows <- nrow(x)
  # can probably show e.g. 60 rows at font size 8
  row.cex <- min(0.5, 60/n_rows)
  
  # draw the heatmap
  g<-Heatmap(x,
          cluster_rows = FALSE,
          col = exprs_cols,
          row_names_gp = gpar(fontsize = row_names_gp,
                              cex = row.cex),
          cluster_columns = FALSE,
          show_column_names = FALSE,
          row_title = NULL,
          row_split=top_markers$cluster,
          column_split = cell_clusters[cell_order],
          top_annotation = subgroupAnnotation,
          right_annotation = clusterAnnotation,
          row_gap=unit(0.5,"mm"),
          column_title_side="bottom",
          column_gap=unit(0.5, "mm"),
          use_raster=TRUE,
          raster_device="CairoPNG",
          raster_quality=4,
          show_heatmap_legend = FALSE
  )
}