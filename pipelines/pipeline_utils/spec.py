import os

def get(infile, outfile, PARAMS):
    '''
    Make a dictionary of commonly need paths and parameteres
    '''

    parts = outfile.split("/")
    print("***********")
    print(parts)
    nparts = len(parts)

    spec = { "sample_name": parts[0].split(".")[0],
             "outdir": os.path.dirname(outfile),
             "indir": os.path.dirname(infile),
             "resolutions" : [x.strip()
                              for x in
                              PARAMS["runspecs_cluster_resolutions"].split(",")]
             }
    if PARAMS["input_format"] == "rds":
        spec["seurat_object"] = os.path.join(parts[0], "begin.rds")
    else:
        spec["seurat_object"] = os.path.join(parts[0], "begin.h5seurat")

    # take care of making the output directory.
    if not os.path.exists(spec["outdir"]):
        os.mkdir(spec["outdir"])

    if outfile.endswith(".sentinel"):
        spec["log_file"] = outfile.replace(".sentinel", ".log")

    # if we are in the component directory
    if nparts >= 3:

        spec["components"] = parts[1].split(".")[1]
        graph_path = os.path.join(parts[0], parts[1], "graph.rds")

        if os.path.exists(graph_path):
            spec["graph"] =  graph_path

    # if we are in the cluster directory
    if nparts >= 4:

        spec["resolution"] = parts[2][len("cluster."):-len(".dir")]

    return(spec)
