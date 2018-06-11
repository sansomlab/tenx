## Title ----
##
## Fetch geneset annotations.
##
## Description ----
##
## This script retrieves ID mappings (Ensembl to entrez)
## and KEGG pathway information.
##
## Details ----
##
## Ensembl to Entrez mappings are retrieved using biomaRt. KEGG pathways are retrieved 
## directly from KEGG.
##
## Usage ----
##
## $ Rscript getGenesetAnnotations.R
##           --ensemblversion=latest
##           --species=mm
##           --outdir=.


# Libraries ----

stopifnot(
  require(optparse),
  require(gsfisher)
)

# Options ----

option_list <- list(
    make_option(
      c("--ensemblversion"), 
      default="latest",
      help="either latest or a specific number"
      ),
    make_option(
      c("--species"), 
      default="none",
      help="species - mm or hs"
      ),
    make_option(
      c("--outdir"), 
      default="none",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

# Fetch Ensembl to Entrez ID mappings ----

anno <- fetchAnnotation(species=opt$species)
write.table(anno,
            gzfile(file.path(opt$outdir,"ensembl.to.entrez.txt.gz")),
            quote=FALSE,
            row.names=FALSE,sep="\t")

# Fetch KEGG pathways ----

kegg_pathways <- fetchKEGG(species=opt$species)
saveRDS(kegg_pathways, file=file.path(opt$outdir,"kegg_pathways.rds"))