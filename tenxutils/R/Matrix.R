#' Slightly modified version of Read10x from Seurat.
#' @param data.dir  The directory containing the 10x matrix, cell names and barcodes.
Read10X <- function(data.dir = NULL){

    if (file.exists(file.path(data.dir,"barcodes.tsv"))) {
        barcodes <- readLines(file.path(data.dir, "barcodes.tsv"))
    } else {
        if(file.exists(file.path(data.dir, "barcodes.tsv.gz"))) {
            barcodes <- readLines(gzfile(file.path(data.dir, "barcodes.tsv.gz")))
        } else {
            stop("Barcode file missing")
        }}

    if(file.exists(file.path(data.dir, "matrix.mtx"))) {
            data <- readMM(file = file.path(data.dir, "matrix.mtx"))
        } else {
            if(file.exists(file.path(data.dir, "matrix.mtx.gz"))) {
            data <- readMM(file = gzfile(file.path(data.dir, "matrix.mtx.gz")))
            } else {
      stop("Matrix file missing")
            }}

     if(file.exists(file.path(data.dir, "features.tsv"))) {
         features = read.table(file.path(data.dir,"features.tsv"), as.is=T)
     } else {
         if(file.exists(file.path(data.dir, "features.tsv.gz")))
         {
             features = read.table(file.path(data.dir, "features.tsv.gz"), as.is=T)
         } else {
             stop("Features file missing")
         }}

    rownames(data)  <- make.unique(
      names = as.character(
        ## TODO: could take the ENSEMBL ID rather than
        ## the symbols
        ## currently we explicity track the new
        ## ID -> symbol mapping (see below)
        x = features$V2
      ))

    colnames(x = data) <- barcodes

    return(data)
}


#' Split barcode and aggregation identifier
#'
#' Cell barcodes follow the format "barcode-aggegationId"
#'
#' @param codes A chracter vector of cell barcodes.
#'
#' @return A data.frame of two columns:
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
#' @param features_tsv_file The path to the original features tsv file.
#' @param metadata Table of cell metadata.
#' @param gzip Should the output files be compresssed?
#' Includes \code{c("code", "agg_id", "barcode", "seq_id", "sample_id"}
#' and any metadata encoded in the sample file names.
#'
#' @return \code{TRUE} if successful.
writeMatrix <- function(
                        dir, matrix, barcodes, features_tsv_file, metadata,
                        gzip = TRUE
){

  if (!dir.exists(dir)) {
    dir.create(dir)
  }

    matrix_path = file.path(dir, "matrix.mtx")
    barcodes_path = file.path(dir, "barcodes.tsv")
    features_path = file.path(dir, basename(features_tsv_file))
    metadata_path = file.path(dir, "metadata.tsv")

    ## write out the matrix
    writeMM(matrix, matrix_path)

    ## write out the "cell" barcodes
    write.table(barcodes, barcodes_path,
                col.names=FALSE, sep=",", row.names=FALSE, quote=FALSE)

    ## copy over the gene names (!)
    file.copy(features_tsv_file, features_path)

    ## write out the metadata
    write.table(
        metadata, metadata_path,
        col.names=TRUE, sep="\t", row.names=FALSE, quote=FALSE)


    if(gzip)
    {
        gzip(matrix_path, overwrite=TRUE)
        gzip(barcodes_path, overwrite=TRUE)
        gzip(metadata_path, overwrite=TRUE)
        if(!endsWith(features_path, ".gz"))
        {
            gzip(features_path)
        }
    } else {
        if(endsWith(features_path, ".gz"))
        {
            gunzip(features_path, overwrite=TRUE)
        }

        }

  # write out the metadata table
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

    allf <- modelGeneVar(sce_object, parametric=TRUE, subset.row = rowMeans(as.matrix(logcounts(sce_object))) > min_mean)

    hvg.out <- allf[which(allf$FDR <= p_adjust_threshold),]

    hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]

    hvg.out
}


#' Run FastExpMean by chunk to conserve memory
#' @param x a matrix
#' @param nrows number of rows per chunk

FastExpMeanChunked <- function(x,
                               rows_per_chunk=2000
                               )
{

  result <- c()

  nrows_x <- nrow(x)

  for(i in seq(1, nrows_x, rows_per_chunk)){

    j <- min(i + rows_per_chunk - 1, nrows_x)

    result <- c(result,
                Seurat:::FastExpMean(Matrix(x[i:j,],sparse=TRUE), display_progress = FALSE))
  }
  result
}



#' Run rowSum by chunk to conserve memory
#' @param x a matrix
#' @param nrows number of rows per chunk

ExpMeanMatrixChunked <- function(x,
                               rows_per_chunk=2000
                               )
{

  result <- c()

  nrows_x <- nrow(x)

  for(i in seq(1, nrows_x, rows_per_chunk)){

    j <- min(i + rows_per_chunk - 1, nrows_x)

    result <- c(result,
                rowSum(Matrix(x[i:j,]))/ncol(x))
  }
  result
}

