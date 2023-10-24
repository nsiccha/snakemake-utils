import json
import sys
import pathlib
import pandas as pd
import logging
from argparse import Namespace


def load_json(path): 
    with open(path, "r") as f: 
        return json.load(f)


def handle_output(path, wildcards_and_params, data):
    path = pathlib.Path(path)
    if path.suffix == ".json":
        assert isinstance(data, dict)
        logging.info(f"Writing to {path}")
        return write_json(path, wildcards_and_params | data)
    if path.suffix == ".csv":
        assert isinstance(data, pd.DataFrame)
        logging.info(f"Writing to {path}")
        return data.to_csv(path)
    raise NotImplementedError(f"Don't know how to handle {path}!")


def write_json(path, data):
    with open(path, 'w') as f: 
        json.dump(data, f)
    return data


def identity(x): return x


def handle_log(log):
    if len(log):
        stdout = log.get("stdout")
        stderr = log.get("stderr")
        if stdout is None and stderr is None:
            stdout = stderr = log[0]
        if stdout == stderr:
            logging.info(f"Redirecting stdout & stderr to {stderr}")
            sys.stdout = sys.stderr = open(stderr, "w")
        else:
            if stdout is not None: 
                logging.info(f"Redirecting stdout to {stdout}")
                sys.stdout = open(stdout, "w")
            if stderr is not None: 
                logging.info(f"Redirecting stdout to {stderr}")
                sys.stderr = open(stderr, "w")


def auto(ctx, fn=None, post=identity):
    if fn is None: 
        return lambda fn: auto(ctx, fn, post=post)

    snakemake = ctx.get("snakemake", Namespace(
        config=dict(),
        log=[],
        wildcards=dict(), params=dict(),
        input=[], output=[]
    ))
    config = snakemake.config
    logging.basicConfig(level=getattr(logging, config.get("logging", "INFO")))
    script_path = ctx.get("__real_file__", ctx.get("__file__"))
    logging.info(f"Running {pathlib.Path(script_path).resolve()}")
    handle_log(snakemake.log)
    wildcards_and_params = dict(**snakemake.wildcards, **snakemake.params)
    logging.info("Calling snakemake function")
    rv = fn(
        input=snakemake.input, 
        output=snakemake.output, 
        **config, **wildcards_and_params
    )
    if rv is not None: 
        logging.info("Automatically writing to outputs")
        rv = [
            handle_output(outputi, wildcards_and_params, rvi)
            for outputi, rvi in zip(snakemake.output, rv)
        ]
    logging.info("Finished successfully")
    return rv