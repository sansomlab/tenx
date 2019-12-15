# Helper functions for working with Seurat objects.

## function to compute the correlation between clusters
## in reduced dimension (e.g. principle component) space
#' @param seurat_object A seurat object containing reduced dimension embeddings and cluster ids.
#' @param white_list A file containing a list of the cell barcodes to retain
#' @param black_list A file containing a list of the cell barcodes to remove
#' @param subset_factor A factor (metadata column) to subset on
#' @param subset_factor_level The level of the subset_factor to be retained
seurat_subset_cells <- function(seurat_object=NULL,
                                 subset_white_list = NULL,
                                 subset_black_list = NULL,
                                 subset_factor = NULL,
                                 subset_factor_level)
{
 
  cells_to_keep <- Cells(s)
  message("original number of cells", length(cells_to_keep))
    
  # 1. subset to whitelist
  if(!is.null(subset_white_list))
  {
    message("subsetting to whitelisted cells")
    whitelist <- scan(subset_white_list, "character")
    cells_to_keep <- intersect(cells_to_keep, whitelist)
    message("No. cells after whitelisting", length(cells_to_keep))
  }
    
  # 2. subset to blacklist
  # Remove blacklisted cells.
  if(!is.null(subset_black_list))
  {
    message("removing blacklisted cells")
    
    blacklist <- scan(subset_black_list, "character")
    cells_to_keep <- cells_to_keep[!cells_to_keep %in% blacklist]
    message("No. cells ids after blacklisting", length(cells_to_keep))
  }
  
  # 3. subset to factor
  if (!is.null(subset_factor)) {
    if(!subset_factor %in% colnames(s[[]])) {
      stop("The given subsetting factor must match a column in the metadata")
    }
    
    if (!subset_factor_level %in% s[[]][,opt$subsetfactor]) {
      stop("The specified level of the subsetting factor does not exist")
    }
    
    message("subsetting by factor")
    cell_subset <- rownames(s[[]])[
      s[[]][, subset_factor] == subset_factor_level
      ]
    cells_to_keep <- intersect(cells_to_keep, cell_subset)
    message("No. cells ids after subsetting to factor", length(cells_to_keep))
  }
    
    message("Number of cells before subsetting:")
    print(length(colnames(x = s)))
    
    s <- SubsetData(s, cells=cells_to_keep)
    
    message("Number of cells after subsetting:")
    print(length(colnames(x = s)))
    
    return(s)
}

## function to compute the correlation between clusters
## in reduced dimension (e.g. principle component) space
#' @param seurat_object A seurat object containing reduced dimension embeddings and cluster ids.
#' @param metadata_file The path to a tsv file containing the metadata
seurat_add_metadata <- function(seurat_object=NULL,
                                metadata_file=NULL)
{
  metadata <- read.table(metadata_file,
                         sep="\t", header=TRUE, as.is=TRUE)

    # ensure that the barcode column is present in the metadata
  if (!"barcode" %in% colnames(metadata)) {
    stop('Mandatory "barcode" column missing from the metadata')
  }
  
  rownames(metadata) <- metadata$barcode
  metadata$barcode <- NULL
  
  metadata <- metadata[colnames(x = s), ]
  
  for(meta_col in colnames(metadata)){
    s[[meta_col]] <- metadata[[meta_col]]
  }

  return(s)
  }

#' Get a breakdown of the number of cells in a seurat object 
#'
#' Initialises a data.frame or appends a new column
#' with a user-defined tag.
#'
#' @param s Seurat object
#' @param cell_numbers data.frame of an earlier call to
#' getCellNumbers, if applicable.
#' @param stage Character value to tag the new entry.
#'
#' @return A data.frame
seurat_track_cell_numbers <- function(s, 
                             cell_numbers=NULL, 
                             stage="input",
                             groupby=opt$groupby) {

  counts <- as.data.frame(table(s[[]][,groupby]))
  colnames(counts) <- c(groupby, stage)
  rownames(counts) <- counts[[groupby]]
  counts[[groupby]] <- NULL
  
  if ( is.null(cell_numbers) ) {
    result <- counts
    colnames(counts) <- stage
  } else {
    cnames <- c(colnames(cell_numbers), stage)
    result <- cbind(cell_numbers, counts[[stage]])
    colnames(result) <- cnames
  }
  return(result)
}



#' Compute percentage of mitochondrial genes and add to a seurat object
#' Compatible with mouse and human data! 
#'
#' Initialises a data.frame or appends a new column
#' with a user-defined tag.
#'
#' @param s Seurat object
seurat_pc_mito <- function(s) {
  ## calculate the percentage of mitochondrial genes here
  ## In humans the pattern is ^MT-, in mouse it is ^mt-: to accomodate both
  ## we simply ignore the case (specificity is still maintained).
  mito.genes <- grep("^MT-", rownames(x = s), value=TRUE, ignore.case=TRUE)

  ## NOTE: unlike in the Seurat vignette, our data is not yet log transformed.
  # TODO: is Matrix:: required here?
  percent.mito <- Matrix::colSums(GetAssayData(object = s)[mito.genes, ]) / Matrix::colSums(GetAssayData(object = s))

  ## AddMetaData adds columns to object@data.info,
  ## and is a great place to stash QC stats
  s$percent.mito <- percent.mito
  
  return(s)
}
  