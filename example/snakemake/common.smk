rule lint:
    input: "snakefile", "snakemake/common.smk"
    output: "output/lint.log"
    shell: "snakemake --lint > {output} 2> {log} || echo 'Check linting output!'"