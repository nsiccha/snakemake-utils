import snakemakeutils

@snakemakeutils.auto(globals())
def test(input, **kwargs):
    yield kwargs   
         