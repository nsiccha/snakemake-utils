import json
import inspect


def load_json(path): 
    with open(path, "r") as f: return json.load(f)

def write_json(path, data):
    with open(path, 'w') as f: json.dump(data, f)
    return data

def identity(x): return x

def auto(fn, snakemake=None, post=identity):
    rv = fn(input=snakemake.input, output=snakemake.output, **snakemake.wildcards)
    if rv is None: return
    return [
        write_json(outputi, post(dict(**snakemake.wildcards) | rvi))
        for outputi, rvi in zip(snakemake.output, rv)
    ]