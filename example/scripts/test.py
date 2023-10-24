import snakemakeutils

print(globals()["snakemake"])

def test(input, **kwargs):
    print(input, kwargs)
    yield kwargs

snakemakeutils.auto(test, snakemake)
