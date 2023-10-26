import snakemakehelpers

@snakemakehelpers.auto(globals())
def test(input, **kwargs):
    if int(kwargs["seed"]) % 2: raise ValueError("I fail if seed % 2!")
    yield kwargs      
         