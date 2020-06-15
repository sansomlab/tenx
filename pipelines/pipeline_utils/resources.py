import math

def get(memory="4G", cpu=1):
    '''calculate the resource requirements and return a
       dictonary that can be used to update the local variables'''

    if not memory.endswith("G"):
        raise ValueError("Memory must be specified as XXG")

    gb_requested = int(memory[:-1])

    mem_gb = int(math.ceil(gb_requested / float(cpu) ))

    spec = {"job_memory": str(mem_gb) + "G",
            "job_threads": cpu,
            "r_memory": gb_requested * 1000}

    return spec
