import pandas
import InquirerPy
import snakemakeutils

@snakemakeutils.auto
def test(input, output, **kwargs):
    print(input, output, kwargs)
    return kwargs