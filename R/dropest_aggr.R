## Aggregate the dropEst matrices

message("cellranger_postprocessAggrMatrix.R")
timestamp()

# Libraries ----

stopifnot(
    require(optparse),
    require(methods), # https://github.com/tudo-r/BatchJobs/issues/27
    require(Matrix),
    require(S4Vectors),
    require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
        c("--sampletable"),
        dest = "sampletable",
        help="Input sample table (sample_id, seq_id, agg_id)"
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


    if(colnames(emat)!=colnames(nmat)) { stop("This is not expected.") }


    all_mat[[sample_id]] <- rbind(emat,nmat)
    colnames(all_mat[[sample_id]]) <- paste(colnames(all_mat[[sample_id]]),
                                            samples[sample_id,"agg_id"],
                                            sep="-")

    emat_names <- unique(c(emat_names, rownames(emat)))
    nmat_names <- unique(c(nmat_names, rownames(nmat)))

    cell_names <- c(cell_names, colnames(all_mat[[sample_id]]))
}

gene_names <- c(emat_names, nmat_names)

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

dir <- opt$outdir

message("writing out the results")
# write out the data matrix
writeMM(results, file.path(dir, "matrix.mtx"))

# write out the "cell" barcodes
write.table(
    data.frame(x=colnames(results)), file.path(dir, "barcodes.tsv"),
    col.names=FALSE, sep=",", row.names=FALSE, quote=FALSE)

# write out the genes
write.table(
    data.frame(x=rownames(results)), file.path(dir, "genes.tsv"),
    col.names=FALSE, sep=",", row.names=FALSE, quote=FALSE)

timestamp()
message("Completed")
