import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvelo as scv
from scipy import io
import os
import logging
import sys
import matplotlib
import argparse

matplotlib.use('Agg')

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("run_scvelo")

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

parser = argparse.ArgumentParser()
parser.add_argument("--dropest_dir", default=1, type=str,
                    help="File with dropest layers")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--cluster_colors", default=1, type=str,
                    help="tsv file with the color palette for the clusters")
parser.add_argument("--rdims", default="1", type=str,
                    help="reduced dimensions file")
parser.add_argument("--rdim_method", default="umap", type=str,
                    help="name of the dimension reduction method")
parser.add_argument("--rdim1", default="UMAP1", type=str,
                    help="reduced dimension 1")
parser.add_argument("--rdim2", default="UMAP2", type=str,
                    help="reduced dimension 2")


args = parser.parse_args()

# ########################################################################### #
# ######################## Initialise AnnData ############################### #
# ########################################################################### #

exon_matrix = os.path.join(args.dropest_dir, "exons.mtx.gz")
intron_matrix = os.path.join(args.dropest_dir, "introns.mtx.gz")
spanning_matrix =  os.path.join(args.dropest_dir, "spanning.mtx.gz")

exons = io.mmread(exon_matrix).transpose().tocsr()
introns = io.mmread(intron_matrix).transpose().tocsr()
spanning = io.mmread(spanning_matrix).transpose().tocsr()

adata = ad.AnnData(X=exons)
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = introns
adata.layers["ambiguous"] = spanning

# Add the variable and observation information
feat = pd.read_csv(os.path.join(args.dropest_dir,"features.tsv.gz"), header=None)
feat.columns = ["gene_symbol"]
feat.index = feat["gene_symbol"]
samples  = pd.read_csv(os.path.join(args.dropest_dir, "barcodes.tsv.gz"),header=None)
metadata = pd.read_csv(os.path.join(args.dropest_dir,"metadata.tsv.gz"), sep="\t")
metadata.index = metadata.barcode.values

adata.vars = feat
adata.obs = metadata

L.info("AnnData initialised")


# ########################################################################### #
# ################### Get rdims (e.g. umap) info  ########################### #
# ########################################################################### #

rdims = pd.read_csv(args.rdims, sep="\t")

rdims.index = rdims["barcode"]

rdata = adata[rdims.index]
rdata.obs["clusters"] = rdims["cluster"].values
rdata.obsm["X_" + args.rdim_method] = rdims[[args.rdim1, args.rdim2]].values

L.info("Rdims info added")


# ########################################################################### #
# ############################ run scvelo ################################### #
# ########################################################################### #

scv.settings.figdir = args.outdir + "/"

# Get the color palette
ggplot_palette = [x for x in pd.read_csv(args.cluster_colors,
                      header=None, sep="\t")[0].values]

# preprocessing
scv.pp.filter_and_normalize(rdata)

scv.pp.moments(rdata)

# velocity computation
scv.tl.velocity(rdata)

scv.tl.velocity_graph(rdata)

L.info("Velocity calculation completed")

# make the plots
scv.pl.velocity_embedding_stream(rdata, basis=args.rdim_method,dpi=300,
                                 palette=ggplot_palette,
                                 save="stream.png")

scv.pl.velocity_embedding(rdata, basis=args.rdim_method,
                          arrow_length=2, arrow_size=1.5,
                          dpi=300,
                          palette=ggplot_palette,
                          save="velocity.png")

scv.pl.velocity_embedding_grid(rdata, basis=args.rdim_method,
                               dpi=300,
                               palette=ggplot_palette,
                               save="grid.png")

L.info("Plotting finished")
