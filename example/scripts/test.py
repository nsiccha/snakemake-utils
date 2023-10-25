import snakemakeutils

@snakemakeutils.auto(globals())
def test(input, **kwargs):
    if kwargs["seed"] % 2: raise ValueError("I fail if seed % 2!")
    yield kwargs      
         