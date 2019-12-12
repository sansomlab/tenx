## 1. Aggregate the dropEst matrices
## 2. Subset the matrices to a given set of barcodes

message("cellranger_postprocessAggrMatrix.R")
timestamp()

# Libraries ----

stopifnot(
    require(optparse),
    require(methods), # https://github.com/tudo-r/BatchJobs/issues/27
    require(Matrix),
    require(S4Vectors),
    require(tenxutils),
    require(R.utils)
)

# Options ----

option_list <- list(
    make_option(
        c("--sampletable"),
        dest = "sampletable",
        help="Input sample table (sample_id, seq_id, agg_id)"
    ),
    make_option(
        c("--barcodes",
          default = NULL,
          dest = "barcodes",
          help = "A tsv file containing the list of barcodes that should be retained")
    ),

    make_option(
        c("--outdir"),
        default=".",
        dest = "outdir",
        help="Location for outputs files. Must exist."
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Functions ----

# Input data ----

samples <- read.table(opt$sampletable, header=T, sep="\t", as.is=T)
rownames(samples) <- samples$sample_id

all_mat <- list()

## unfotunately rownames are not consistent so we have to operate
## over the data twice

emat_names <- c()
nmat_names <- c()
smat_names <- c()
cell_names <- c()
for(sample_id in samples$sample_id)
{
    m_file <- file.path(paste0(sample_id, "-dropest"),
                        paste0(sample_id, ".matrices.rds"))

    m <- readRDS(m_file)

    emat <- m$exon
    rownames(emat) <- paste0("ex_", rownames(emat))

    nmat <- m$intron
    rownames(nmat) <- paste0("in_", rownames(nmat))

    smat <- m$spanning
    rownames(smat) <- paste0("sp_", rownames(smat))

    if(colnames(emat)!=colnames(nmat)) { stop("This is not expected.") }


    all_mat[[sample_id]] <- rbind(emat,nmat,smat)

    colnames(all_mat[[sample_id]]) <- paste(colnames(all_mat[[sample_id]]),
                                            samples[sample_id,"agg_id"],
                                            sep="-")

    emat_names <- unique(c(emat_names, rownames(emat)))
    nmat_names <- unique(c(nmat_names, rownames(nmat)))
    smat_names <- unique(c(smat_names, rownames(smat)))

    cell_names <- c(cell_names, colnames(all_mat[[sample_id]]))
}

gene_names <- c(emat_names, nmat_names, smat_names)

results <- Matrix(data=0,
                  length(gene_names),
                  length(cell_names))

rownames(results) <- gene_names
colnames(results) <- cell_names

message("building results matrix")
for(sample in names(all_mat))
{
    results[rownames(all_mat[[sample]]),
            colnames(all_mat[[sample]])] <- all_mat[[sample]]
}


## Perform subsetting
if(!is.null(opt$barcodes)) {

    message("Subsetting to given list of barcodes")

    ## subset to the given list of barcodes
    barcodes_to_keep <- read.table(gzfile(opt$barcodes))$V1

    message("Length of the given list of barcodes: ", length(barcodes_to_keep))

    message("number of barcodes before subsetting: ", ncol(results))
    results <- results[,colnames(results) %in% barcodes_to_keep]

    message("number of barcodes after subsetting: ", ncol(results))
}


print(dim(results))

message("writing out the results")
matrixFile <- file.path(opt$outdir, "matrix.mtx")

## write out the data matrix
## (writeMM does not support writing to a connection)

writeMM(results, matrixFile)
gzip(matrixFile, overwrite = TRUE)

## write out the "cell" barcodes
barcodeFile <- gzfile(file.path(opt$outdir, "barcodes.tsv.gz"), "wt")
write.table(
    data.frame(x=colnames(results)), barcodeFile,
    col.names=FALSE, sep=",", row.names=FALSE, quote=FALSE)
close(barcodeFile)

## write out the features
featureFile <- gzfile(file.path(opt$outdir, "features.tsv.gz"), "wt")
write.table(
    data.frame(x=rownames(results)), featureFile,
    col.names=FALSE, sep=",", row.names=FALSE, quote=FALSE)
close(featureFile)

timestamp()
message("Completed")
