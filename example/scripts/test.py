import snakemakeutils
import time

print(__real_file__)
print(globals()["snakemake"])
print(dict(snakemake.wildcards) | snakemake.config)

print(dir(snakemake))

for name in ['bench_iteration', 'config', 'input', 'log', 'log_fmt_shell', 'output', 'params', 'resources', 'rule', 'scriptdir', 'threads', 'wildcards']:
    attr = getattr(snakemake, name)
    print(name, type(attr), attr)
    try:
        print(list(attr))
        print(dict(attr))
    except:
        print("xxxx")
    print("========")
# print(dict(snakemake.rule))
# print(dict(snakemake.log))
# print(dict(snakemake.log))



# @snakemakeutils.auto(globals())
# def test(input, ctx=globals(), **kwargs):
#     print(ctx)
#     print(input, kwargs)
#     yield kwargs   
