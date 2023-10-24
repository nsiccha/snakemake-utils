import inspect

def identity(x): return x

def auto_snakemake(snakemake, fn, post=identity):
    rv = fn(snakemake.input, snakemake.output, **snakemake.wildcards)
    if rv is None: return
    write_json(snakemake.output[0], post(dict(**snakemake.wildcards) | rv))
    return rv