import os
import re
import argparse
import anndata
import pandas as pd
import logging
import sys



# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)


# ########################################################################### #
# ######################## Parse the arguments ############################## #
# ########################################################################### #

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--reduced_dims_matrix_file", default="reduced_dims.tsv.gz", type=str,
                    help="File with reduced dimensions")
parser.add_argument("--barcode_file", default="barcodes.tsv.gz", type=str,
                    help="File with the cell barcodes")
parser.add_argument("--outdir",default=1, type=str,
                    help="path to output directory")
parser.add_argument("--comps", default="1", type=str,
                    help="Number of dimensions to include in the knn computation")
parser.add_argument("--method", default="scanpy", type=str,
                    help="scanpy|hnsw|precomputed (scanpy uses pynndescent)")
parser.add_argument("--knn", default="none", type=str,
		    help="path to anndata with precomputed knn")
parser.add_argument("--k", default=20, type=int,
                    help="number of neighbors")
parser.add_argument("--metric", default="euclidean", type=str,
                    help="the distance metric")
parser.add_argument("--threads", default=4, type=int,
                    help="number of threads")
parser.add_argument("--fullspeed", default=False, action="store_true",
                    help="number of threads")

args = parser.parse_args()

L.info("Running with arguments:")
print(args)


# ########################################################################### #
# ############## Create outdir and set results file ######################### #
# ########################################################################### #

# write folder
results_file = args.outdir + "/" + "paga_anndata.h5ad"


# ########################################################################### #
# ######################### Make the anndata object ######################### #
# ########################################################################### #

# Read matrix of reduced dimensions, create anndata and add dimensions
reduced_dims_mat = pd.read_csv(args.reduced_dims_matrix_file, sep="\t")
reduced_dims_mat.index = [x for x in pd.read_csv(args.barcode_file, header=None)[0]]

adata = anndata.AnnData(obs=[x for x in reduced_dims_mat.index])
adata.obs.index = reduced_dims_mat.index
adata.obs.rename(columns={0:'barcode'}, inplace=True)

# Add dimensions to anndata
colnames = list(reduced_dims_mat.columns)
colname_prefix = list(set([re.sub(r'_[0-9]+', '', i) for i in colnames]))[0] + "_"
select_comps = args.comps.split(",")
select_comps = [ colname_prefix + item for item in select_comps]
reduced_dims_mat = reduced_dims_mat[select_comps]

n_pcs = len(select_comps)

L.info("Using comps " + ', '.join(list(reduced_dims_mat.columns)))

adata.obsm['X_rdims'] = reduced_dims_mat.to_numpy(dtype="float32")


# ########################################################################### #
# ####################### Nearest neighbor computation ###################### #
# ########################################################################### #

# Run neighbors
L.info( "Using " + str(args.k) + " neighbors")

if args.method == "scanpy":

    L.info("Computing neighbors using default scanpy method")
    from scanpy.preprocessing import neighbors

    neighbors(adata,
              n_neighbors = args.k,
              metric = args.metric,
              use_rep = 'X_rdims')


elif args.method == "hnsw":

    L.info("Computing neighbors using hnswlib (with scvelo a la pegasus!)")
    # we use the neighbors function from scvelo (thanks!)
    # with parameters from pegasus (for a more exact result).

    from scvelo.pp import neighbors

    num_threads = (args.threads if args.fullspeed else 1)

    neighbors(adata,
              n_neighbors = args.k,
              n_pcs = None,
              use_rep = "X_rdims",
              knn = True,
              random_state = 0,
              method = 'hnsw',
              metric = args.metric,
              metric_kwds = {"M":20,
                             "ef":200,
                             "ef_construction":200},
              num_threads=num_threads)

elif args.method == "precomputed":

    L.info("Reading neighbors from precomputed anndata object")  
    #read anndata object with knn graph
    adata_path = args.knn
    L.info("Anndata object:"+ adata_path)
    adata_knn = anndata.read_h5ad(adata_path)
    #add knn to adata
    adata.uns["neighbors"] = adata_knn.uns["neighbors"]
    adata.obsp = adata_knn.obsp
    #rename rep_use to 'X_rdims'
    from scanpy._utils import NeighborsView
    neighbors = NeighborsView(adata, "neighbors")

    print(neighbors["params"])

    neighbors["params"]["use_rep"]="X_rdims"

    adata.uns["neighbors"]["params"]=neighbors["params"]


else:
    raise ValueError("nn method not recognised")

# compute clusters
adata.write(os.path.join(args.outdir,
                         "anndata.h5ad"))

L.info("Complete")
