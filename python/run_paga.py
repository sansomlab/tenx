import os
import argparse
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as pl
import scanpy as sc
import pandas as pd
from scipy import sparse
import pyreadr
#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()


###### Parse arguments ######
#############################

parser = argparse.ArgumentParser()
parser.add_argument("--pcs", default=1, type=str, help="File with pcs")
parser.add_argument("--outdir", default=1, type=str, help="path to output directory")
parser.add_argument("--cluster_ids", default=1, type=str, help="RDS files with cluster ids")
parser.add_argument("--comps", default="1", type=str, help="Number of PCs to include in knn and umap computation")
parser.add_argument("--resolution", default=1, type=str, help="cluster resolution")

args = parser.parse_args()
pcs = args.pcs
outdir = args.outdir
comps = args.comps
resolution = args.resolution
clusterids = args.cluster_ids


###### Create outdir and set results file ######
################################################

# write folder
results_file = outdir + "/" + "paga_anndata.h5ad"

# figures folder
sc.settings.figdir = outdir

###### Run paga ######
######################

# Read pcs, create anndata and add pcs
pcs = pd.read_table(pcs)
adata = sc.AnnData(obs=pcs.index)
adata.obs.index = pcs.index
adata.obs.rename(columns={0:'barcode'}, inplace=True)

# Add pcs
select_comps = comps.split(",")
select_comps = [ "PC_" + item for item in select_comps]
pcs = pcs[select_comps]
message = "Using comps " + ', '.join(list(pcs.columns))
print(message)
adata.obsm['X_pca'] = pcs.to_numpy(dtype="float32")

# Read and add cluster ids
df = pyreadr.read_r(clusterids)
df = df[None][['cluster_id']]
adata.obs['cluster_id'] = df.values

# Run neighbors
sc.pp.neighbors(adata, n_neighbors=20, use_rep='X_pca')

# Run, plot and save umap
sc.tl.umap(adata)
sc.pl.umap(adata, color="cluster_id", legend_loc='on data', save = True, show=False)
adata.obsm["X_umap_init_pre"] = adata.obsm["X_umap"]

# Run and plot paga
sc.tl.paga(adata, groups='cluster_id')
sc.pl.paga(adata, save=True, show=False)

# Run, plot and store paga-initialised umap
sc.tl.umap(adata, init_pos = sc.tl._utils.get_init_pos_from_paga(adata))
sc.pl.umap(adata, color="cluster_id", legend_loc='on data', save = "_initipos", show=False)
adata.obsm['X_umap_init_pos'] = adata.obsm['X_umap']

# Rename X_umap_init_pos
del adata.obsm['X_umap']
adata.obsm['X_umap'] = adata.obsm['X_umap_init_pre']
del adata.obsm['X_umap_init_pre']

# Save paga-initialised UMAP coordinates
umap_pos=pd.DataFrame(adata.obs['barcode'])
umap_pos['x coordinate'] = pd.DataFrame(adata.obsm['X_umap_init_pos'][:,0]).values
umap_pos['y coordinate'] = pd.DataFrame(adata.obsm['X_umap_init_pos'][:,1]).values
out = outdir + '/UMAP_init_pos.csv'
umap_pos.to_csv(out,index=False)

# Write output file
adata.write(results_file)




