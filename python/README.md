# snakemakehelpers

Meant to work with [poetry](https://python-poetry.org/) and [snakemake](https://snakemake.readthedocs.io/en/stable/).

## Installation

Add to an initialized (via `poetry init`) poetry project from PyPI via `poetry add snakemakehelpers`. 

## usage

Run `poetry run snakemake-interact` from a directory with a snakefile to interact with snakemake.

Current options to interact are:

* Process latest snakemake log file:

    Searches the `.snakemake/log` directory for log files and parses the most recent one.
* Select target(s) to make (poetry run snakemake -Fn):

    Generates a list of possible output files via `poetry run snakemake -Fn`, from which the user can choose one or more to make.
* Make target(s) ( poetry run snakemake ... -c):

    Makes user provided targets. 
* Inspect selected snakemake log file (...):

    Inspects the selected log file using `less`. 
* Select a different snakemake log file:

    Select another log file from `.snakemake/log` to process.
* Print internal state of this program