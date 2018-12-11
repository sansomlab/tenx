#' Slightly modified version of Read10x from Seurat.
#' @param data.dir  The directory containing the 10x matrix, cell names and barcodes.
Read10X <- function(data.dir = NULL){
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (! dir.exists(run)){
      stop("Directory provided does not exist")
    }
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "genes.tsv")
    matrix.loc <- paste0(run, "matrix.mtx")
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)

    rownames(x = data) <- make.unique(
      names = as.character(
        ## TODO: could take the ENSEMBL ID rather than
        ## the symbols
        ## currently we explicity track the new
        ## ID -> symbol mapping (see below)
        x = read.table(gene.loc)$V2
      ))

    if (is.null(x = names(x = data.dir))) {
      if(i < 2){
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}


#' Split barcode and aggregation identifier
#'
#' Cell barcodes follow the format "barcode-aggegationId"
#'
#' @param codes A chracter vector of cell barcodes.
#'
#' @return A data.frame of two columns:
#' \described{
#'   \item{code}{The sequenced cell barcode.}
#'   \item{agg_id}{The sample aggregation identifier.}
#' }
barcode2table <- function(codes){
  barcodeTable <- read.table(text=codes, sep="-", as.is = TRUE)
  colnames(barcodeTable) <- c("code","agg_id")
  rownames(barcodeTable) <- codes
  barcodeTable$agg_id <- as.factor(barcodeTable$agg_id)
  return(barcodeTable)
}


#' Write count data in the 10x format, with additional metadata
#'
#' Create a directory containing the count matrix and cell/gene annotation
#' from a sparse matrix of UMI counts,
#' in the format produced by the CellRanger software suite.
#' In addition, write an accompanying metadata file.
#'
#' @param dir Path to the output directory.
#' (Does not need to exist yet).
#' @param matrix A sparse numeric matrix of UMI counts.
#' @param barcodes A character vector of cell barcodes,
#' one per column of \code{x}.
#' @param gene_tsv_file The path to the original genes.tsv file.
#' @param metadata Table of cell metadata.
#' Includes \code{c("code", "agg_id", "barcode", "seq_id", "sample_id"}
#' and any metadata encoded in the sample file names.
#'
#' @return \code{TRUE} if successful.
writeMatrix <- function(
  dir, matrix, barcodes, gene_tsv_file, metadata
){

  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  # write out the data matrix
  writeMM(matrix, file.path(dir, "matrix.mtx"))

  # write out the "cell" barcodes
  write.table(
    barcodes, file.path(dir, "barcodes.tsv"),
    col.names=FALSE, sep=",", row.names=FALSE, quote=FALSE)

  # copy over the gene names (!)
  file.copy(gene_tsv_file, file.path(dir, "genes.tsv"))

  # write out the metadata table
  write.table(
    metadata, file.path(dir, "metadata.tsv"),
    col.names=TRUE, sep="\t", row.names=FALSE, quote=FALSE)

  return(TRUE)
}

#' Downsample a matrix of UMI counts,
#' @param matrixUMI sparse matrix of umi counts
#' @param downsample_method the statistic to normalise per-library
#' @param library_ids a vector of the library_ids
downsampleMatrix <- function(matrixUMI, downsample_method="median", library_ids=c())
{
  requireNamespace("DropletUtils")
  # downsampleCounts drops the colnames
  backup_colnames <- colnames(matrixUMI)

  # Total UMI count per cell
  nUMIs <- Matrix::colSums(matrixUMI)

  # Summary UMI metric by sample (e.g. median)

  statUMIs <- tapply(nUMIs, library_ids, downsample_method)
  cat(sprintf("%s(UMIs) before downsampling (by agg_id):\n", downsample_method))
  print(statUMIs)

  # Identify proportion to retain for each sample
  normFactorUMIs <- min(statUMIs) / statUMIs
  cat("Proportion for downsampling by agg_id:\n")
  print(normFactorUMIs)

  # Extend scaling factor to all cells by sample
  colNormFactor <- normFactorUMIs[library_ids]

  # Downsample
  matrixUMI <- DropletUtils::downsampleMatrix(as(matrixUMI, "dgCMatrix"), colNormFactor, bycol=TRUE)

  # downsampleCounts drops the colnames
  colnames(matrixUMI) <- backup_colnames

  # Ensure the downsampling has worked properly
  nUMIs <- Matrix::colSums(matrixUMI)
  statUMIs <- tapply(nUMIs, library_ids, downsample_method)
  cat(sprintf("%s(UMIs) after downsampling (by agg_id):\n", downsample_method))
  print(statUMIs)
  matrixUMI
}

#' Get highly variable genes using the trend var method from scran
#' @param seurat_object a seurat object
#' @param min_mean minimum mean counts threshold
#' @param p_adjust_threshold threshold for significance
getHVG <- function(seurat_object,
                     min_mean=0,
                     p_adjust_threshold=0.05)
{
    require(scran)
    require(scater)
    sce_object <- as.SingleCellExperiment(seurat_object)

    var.fit.nospike <- trendVar(sce_object,
                                parametric=TRUE,
                                use.spikes=FALSE,
                                loess.args=list(span=0.2))

    var.out.nospike <- decomposeVar(sce_object,
                                    var.fit.nospike,
                                    subset.row=rowMeans(
                                        as.matrix(logcounts(sce_object))) > min_mean)

    ## TODO: some plots should be made.
    ## plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6,
    ##     xlab="Mean log-expression", ylab="Variance of log-expression")
    ## curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
    ## points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)
    hvg.out <- var.out.nospike[which(var.out.nospike$FDR <= p_adjust_threshold),]

    hvg.out <- hvg.out[order(hvg.out$bio,
                             decreasing=TRUE),]

    hvg.out
}
