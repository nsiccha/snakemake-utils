import snakemakeutils

print(globals()["snakemake"])

def test(input, output, **kwargs):
    print(input, output, kwargs)
    yield kwargs

snakemakeutils.auto(test, snakemake)
