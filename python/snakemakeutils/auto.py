import json
import sys
import pathlib

def load_json(path): 
    with open(path, "r") as f: return json.load(f)

def handle_output(path, data):
    path = pathlib.Path(path)
    if path.suffix == ".json":
        return write_json(path, data)
    if path.suffix == ".csv":
        return write_csv(path, data)

def write_json(path, data):
    with open(path, 'w') as f: json.dump(data, f)
    return data

def identity(x): return x

def handle_log(log):
    if len(log):
        stdout = log.get("stdout")
        stderr = log.get("stderr")
        if stdout is None and stderr is None:
            stdout = stderr = log[0]
        if stdout == stderr:
            sys.stdout = sys.stderr = open(stderr, "w")
        else:
            if stdout is not None: 
                sys.stdout = open(stdout, "w")
            if stderr is not None: sys.stderr = open(stderr, "w")

def auto(ctx, fn=None, post=identity):
    if fn is None: return lambda fn: auto(ctx, fn, post=post)

    snakemake = ctx["snakemake"]
    handle_log(snakemake.log)
    rv = fn(input=snakemake.input, output=snakemake.output, **snakemake.config, **snakemake.wildcards, **snakemake.params)
    if rv is None: return
    return [
        write_json(outputi, post(dict(**snakemake.wildcards) | rvi))
        for outputi, rvi in zip(snakemake.output, rv)
    ]