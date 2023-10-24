import snakemakeutils
import time

print(__real_file__)
print(globals()["snakemake"])
print(dict(snakemake.wildcards) | snakemake.config)

@snakemakeutils.auto(globals())
def test(input, ctx=globals(), **kwargs):
    print(ctx)
    print(input, kwargs)
    yield kwargs   
