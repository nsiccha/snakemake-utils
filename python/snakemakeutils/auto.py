import json
import inspect


def load_json(path): 
    with open(path, "r") as f: return json.load(f)

def write_json(path, data):
    with open(path, 'w') as f: json.dump(data, f)

def identity(x): return x

def auto(fn, snakemake=None, post=identity):
    if snakemake is None: snakemake = globals()["snakemake"]
    rv = fn(snakemake.input, snakemake.output, **snakemake.wildcards)
    if rv is None: return
    write_json(snakemake.output[0], post(dict(**snakemake.wildcards) | rv))
    return rv