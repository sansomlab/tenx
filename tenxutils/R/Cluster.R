## function to compute the correlation between clusters
## in principle component space
#' @param seurat_object A seurat object containing pca embeddings and cluster ids.
#' @param pcs The principle components to use.
#' @param cor_method The correlation method to use.
#' @param cluster_average If false, the median pairwise correlation is used. 
clusterCor <- function(seurat_object=NULL,
                       pcs=NULL,
                       cor_method="pearson",
                       cluster_average=FALSE)
{
  pcomps <- t(s@dr$pca@cell.embeddings)[pcs,]
  clusters = as.numeric(as.vector(unique(s@ident)))
  clusters <- clusters[order(clusters)]
  names(clusters) <- paste0("C",clusters)
  n=length(clusters)

  if(cluster_average==FALSE)
  {
    # compute pairwise correlations
    rmat <- matrix(ncol=n, nrow=n)

    for(i in 1:n)
    {
      xclust = clusters[i]
      x <- pcomps[,names(s@ident)[s@ident==xclust]]
      for(j in i:n)
      {
        yclust = clusters[j]
        y <- pcomps[,names(s@ident)[s@ident==yclust]]
        pairwise_cors <- cor(x,y, method=cor_method)
        med_cor <-  median(pairwise_cors)
        rmat[i,j] <- med_cor
        rmat[j,i] <- med_cor
      }
    }
  } else {
    
   # compute correlation of cluster average 
    res <- data.frame(row.names = rownames(pcomps))
    
    # get the cluster averages
    for(i in 1:n)
    {
      xclust <- clusters[i]
      xname <- names(clusters)[i]
      x <- apply(pcomps[,names(s@ident)[s@ident==xclust]],1,mean)
      res[[xname]] <- x[rownames(res)]
    }
    
    rmat <- cor(res, method=cor_method)
    
  }
  rownames(rmat) <- names(clusters)
  colnames(rmat) <- names(clusters)
  rmat
}

