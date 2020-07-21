## Modified FindVariableFeatures to include trendvar
#' @param seurat_object A seurat object
#' @param method One of "trendvar" or other supported by FindVariableFeatures
#' @param nfeatures Number of variable features to find (Seurat methods)
#' @param xlowcutoff For Seurat methods
#' @param xhigcutoff For Seurat methods
#' @param sdcutoff For Seurat methods
#' @param minmean For trendvar
#' @param padj adjusted P value cutoff for trendvar
FindVariableFeaturesMod <- function(seurat_object,
                                    method,
                                    nfeatures,
                                    xlowcutoff,
                                    xhighcutoff,
                                    sdcutoff,
                                    minmean,
                                    padj)
{

    message("Finding variable features")
    ## We need to run FindVariableFeatures to set HVFInfo(object = s)
    ## even if "trendvar" method is specified...
    if(method=="trendvar")
    {
        fvg_method="mean.var.plot"
    } else {
        fvg_method=method
    }

    seurat_object <- FindVariableFeatures(seurat_object,
                                          selection.method=fvg_method,
                                          nfeatures=nfeatures,
                                          mean.cutoff = c(xlowcutoff,
                                                          xhighcutoff),
                                          dispersion.cutoff=c(sdcutoff, Inf))

    xthreshold <- xlowcutoff

    if(method == "trendvar")
    {
        message("setting variable genes using the trendvar method")

        ## get highly variable genes using the getHVG function in tenxutils (Matrix.R)
        ## that wraps the trendVar method from scran.
        hvg.out <- getHVG(seurat_object,
                          min_mean=minmean,
                          p_adjust_threshold=padj)

        ## overwrite the slot
        VariableFeatures(object = seurat_object) <- row.names(hvg.out)

        xthreshold <- minmean
    }
    seurat_object
}

## A function for converting a list of human gene symbols to mouse
#'
#' see: https://github.com/satijalab/seurat/issues/2493
#' and: https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
#'
#' @param human_genes A vector of human gene symbols
convertHumanGeneList <- function(human_gene_symbols)
{
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    genesV2 = getLDS(attributes = c("hgnc_symbol"),
                     filters = "hgnc_symbol",
                     values = human_gene_symbols , mart = human,
                     attributesL = c("mgi_symbol"),
                     martL = mouse,
                     uniqueRows=T)

    mouse_gene_symbols <- unique(genesV2[, 2])

    ## Print the first 6 genes found to the screen
    print(head(mouse_gene_symbols))
    return(mouse_gene_symbols)
}

## Get variable genes from files or from Seurat
#' If files are passed both s and g2m phase lists must be supplied.
#'  - text files with a single headerless column of symbols are
#'    expected.
#' If files are not passed, the lists from Seurat will be used.
#' If species is "mm" and lists from Seurat are being used, the
#' symbols will be converted with bioMart.
#'
#' @param sgenes_file A text file with a single column of gene symbols
#' @param g2mgenes_file A text file with a single column of gene symbols
#' @param species Either "hs" or "mm" - ignored if files are specified
getCellCycleGenes <- function(sgenes_file = NULL,
                              g2mgenes_file = NULL,
                              species = c("mm","hs"))
{
    if(!is.null(sgenes_file) | !is.null(g2mgenes_file))
    {
        message("fetching cell cycle gene symbols from files")
        if(!file.exists(sgenes_file)) {
            stop("sgenes_file does not exist")
        }
        sgenes <- read.table(sgenes_file, header=F, as.is=T)$V1
        if(!file.exists(g2mgenes_file)) {
            stop("g2mgenes_file does not exist")
        }
        g2mgenes <- read.table(g2mgenes_file, header=F, as.is=T)$V1
    } else {

        require(Seurat)

        if(!opt$species %in% c("mm", "hs"))
        { stop ("species not recognised") }

        if(opt$species == "hs") {
            sgenes <- cc.genes$s.genes
            g2mgenes <- cc.genes$g2m.genes
        } else if (opt$species == "mm")
        {
            sgenes <- convertHumanGeneList(cc.genes$s.genes)
            g2mgenes <- convertHumanGeneList(cc.genes$g2m.genes)

        }
    }
    list(s.genes = sgenes,
         g2m.genes = g2mgenes)
}


## Function to load seurat object from either rds or h5seurat
#' This function can load the Seurat object from either rds or 
#' h5seurat format 
#'
#' @param path The path to the stored object
#' @param format A string indicating whether it is 'rds' or 'h5seurat'. If left empty (defaults to "none"), the ending of the path is used to determine file format.
loadSeurat <- function(path,
                       format = "none")
{
    if (endsWith(path, ".rds") | format == "rds") {
        message(sprintf("readRDS: %s", opt$seuratobject))
        s <- readRDS(path)
    } else if (endsWith(path, ".h5seurat") | format == "h5seurat") {
        message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
        stopifnot(require(SeuratDisk))
        s <- LoadH5Seurat(path)
    } else {
        stop("Input format not supported. The format or path ending needs to be rds/.rds or h5seurat/.h5seurat.")
    }
}

## Function to save seurat objects to either rds or h5seurat
#' This function can save the Seurat object to either rds or 
#' h5seurat format 
#'
#' @param path The path where the object should be stored
#' @param format A string indicating whether it is 'rds' or 'h5seurat'. If left empty (defaults to "none"), the ending of the path is used to determine file format.

saveSeurat <- function(path, 
                       format = "none")
{
    if (endsWith(path, ".rds") | format == "rds") {
        message(sprintf("readRDS: %s", opt$seuratobject))
        saveRDS(s, file=file.path(opt$outdir, "begin.rds"))
    } else if (endsWith(path, ".h5seurat") | format == "h5seurat") {
        message(sprintf("LoadH5Seurat: %s", opt$seuratobject))
        stopifnot(require(SeuratDisk))
        SaveH5Seurat(s, file=file.path(opt$outdir, "begin.h5seurat"), overwrite = TRUE)
    } else {
        stop("Output format not supported. The format or path ending needs to be rds/.rds or h5seurat/.h5seurat.") 
    }
}



