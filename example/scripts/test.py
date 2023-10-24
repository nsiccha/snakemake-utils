import snakemakeutils

print(globals()["snakemake"])
print(dict(snakemake.wildcards) | snakemake.config)


def test(input, **kwargs):
    print(input, kwargs)
    yield kwargs

snakemakeutils.auto(test, snakemake)
