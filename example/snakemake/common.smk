rule lint:
    input: "snakefile", "snakemake/common.smk"
    output: "output/lint.log"
    log: "output/lint.log"
    # conda: "env.yaml"
    shell: "snakemake --lint > {output} 2> {log} || echo 'Check linting output!'"