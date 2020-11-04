import os
import math
import types
import pandas as pd

def get_resources(memory="4G", cpu=1):
    '''calculate the resource requirements and return a
       dictonary that can be used to update the local variables'''

    if not memory.endswith("G"):
        raise ValueError("Memory must be specified as XXG")

    gb_requested = int(memory[:-1])

    mem_gb = int(math.ceil(gb_requested / float(cpu) ))

    spec = {"job_memory": str(mem_gb) + "G",
            "job_threads": cpu,
            "r_memory": gb_requested * 1000}

    return(spec["job_threads"],
           spec["job_memory"],
           spec["r_memory"])


def get_vars(infile, outfile, PARAMS, make_outdir=True):
    '''
    Make a dictionary of commonly need paths and parameteres
    '''

    outfile = os.path.relpath(outfile)

    parts = outfile.split("/")

    nparts = len(parts)

    SPEC = { "sample_name": parts[0].split(".")[0],
             "sample_dir": parts[0],
             "seurat_object": os.path.join(parts[0], "begin.rds"),
             "outdir": os.path.dirname(outfile),
             "outname": os.path.basename(outfile),
             }


    if PARAMS["runspecs_cluster_resolutions"]:

        SPEC["resolutions"] = [x.strip()
                               for x in
                               str(PARAMS["runspecs_cluster_resolutions"]).split(",")]

    else:
        SPEC["resolutions"] = []

    if PARAMS["runspecs_predefined_clusters"]:
        SPEC["resolutions"].append("predefined")


    if not infile is None:
        infile = os.path.relpath(infile)
        SPEC["indir"] = os.path.dirname(infile)
        SPEC["inname"] = os.path.basename(infile)

    if(make_outdir):

        # take care of making the output directory.
        if not os.path.exists(SPEC["outdir"]):
            os.makedirs(SPEC["outdir"])

    if outfile.endswith(".sentinel"):
        SPEC["log_file"] = outfile.replace(".sentinel", ".log")

    # if we are in the component directory
    if nparts >= 3 and parts[1].startswith("components."):

        SPEC["component_dir"] = os.path.join(parts[0], parts[1])
        SPEC["components"] = parts[1].split(".")[1]
        graph_path = os.path.join(parts[0], parts[1], "graphs.rds")

        SPEC["seurat_graphs"] =  graph_path

        anndata_path = os.path.join(SPEC["component_dir"],
                                    "anndata.dir",
                                    "anndata.h5ad")

        if os.path.exists(anndata_path):
            SPEC["anndata"] = anndata_path

    # if we are in the cluster directory
    if nparts >= 4 and parts[2].startswith("cluster."):

        SPEC["cluster_dir"] = os.path.join(parts[0], parts[1], parts[2])
        SPEC["resolution"] = parts[2][len("cluster."):-len(".dir")]
        SPEC["cluster_ids"] = os.path.join(SPEC["cluster_dir"], "cluster_ids.rds")
        SPEC["cluster_assignments"] = os.path.join(SPEC["cluster_dir"], "cluster_assignments.tsv.gz")
        SPEC["cluster_colors"] = os.path.join(SPEC["cluster_dir"], "cluster_colors.tsv")

        # get the list of clusters
        cluster_table = os.path.join(SPEC["cluster_dir"], "cluster_ids.tsv")
        if os.path.exists(cluster_table):
            clusters = pd.read_table(cluster_table, header=None)
            clusters = clusters[clusters.columns[0]].tolist()
            SPEC["cluster_table"] = cluster_table
            SPEC["clusters"] = [x for x in set(clusters)]
            SPEC["nclusters"] = len(SPEC["clusters"])

    return([types.SimpleNamespace(**SPEC), SPEC])
