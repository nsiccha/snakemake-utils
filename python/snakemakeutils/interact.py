import pathlib, subprocess, re, io
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
import pandas as pd
import os, shlex
import argparse
from argparse import Namespace
import socket

def print_and_run(cmd):
    print(cmd)
    os.system(cmd)

def mtime(path): return pathlib.Path(path).stat().st_mtime

def jobid(path): return int(re.findall(r"([0-9]+)\.log", path)[0])

def path_to_dict(path): return dict(PATH=path, JOBID=jobid(path), STATE="COMPLETED?") | {
    key: "" for key in ["PARTITION", "NAME", "TIME", "START_TIME", "_7"]
}

def update_logs(info):
    info.log_paths = []
    info.log_path = None
    if not info.path.exists(): 
        return print(f"{info.path} does not exist!")
    base = info.path / ".snakemake"
    if not base.exists():
        return print(f"Snakemake path {base} does not exist!")
    logs_dir = base / "log"
    info.log_paths = sorted((logs_dir).glob("*.log"))
    if not len(info.log_paths):
        return print(f"No snakemake logs exist in {logs_dir}!")
    info.log_path = info.log_paths[-1]
    process_log(info)
    

def process_log(info):
    with open(info.log_path, "r") as fd:
        info.log_content = fd.read()
    info.job_stats = re.search(
        r"Job stats:([\s\S]+?)Select jobs to execute...",
        info.log_content
    )
    if info.job_stats is not None:
        info.job_stats = info.job_stats.group(1).strip()
        # info.job_stats = pd.read_csv(
        #     io.StringIO(info.job_stats.group(1)), sep="\s+"
        # )

    info.output_files = re.findall(
        r"    output: (.+)", 
        info.log_content
    )
    info.slurm_log_paths = re.findall(
        r"Job .+ has been submitted with SLURM jobid .+ \(log: (.+)\).", 
        info.log_content
    )
    info.error_log_paths = re.findall(
        r"log: (.+) \(check log file\(s\) for error details\)", 
        info.log_content
    )
    info.progress = (
        [None] + re.findall(r"[0-9]+ of [0-9]+ steps \(.+\) done", info.log_content)
    )[-1]
    info.jobs_df = None
    if not info.slurm_log_paths: 
        return print("[Jobs] No slurm logs found!")
    info.jobs_df = pd.DataFrame(map(path_to_dict, info.slurm_log_paths))
    for error_log_path in info.error_log_paths:
        info.jobs_df.loc[info.jobs_df.PATH == error_log_path, "STATE"] = "ERRORED"
    info.slurm_queue = subprocess.run(['slurm', 'q'], stdout=subprocess.PIPE).stdout.decode('utf-8')
    for row in pd.read_csv(io.StringIO(info.slurm_queue), sep="\s+").itertuples():
        for key, value in row._asdict().items():
            if key == "Index": continue
            info.jobs_df.loc[info.jobs_df.JOBID == row.JOBID, key] = value
    print(info.jobs_df)

    
def inspect_log(info):
    os.system(f"less {shlex.quote(str(info.log_path))}")

def select_log(info):
    update_logs(info)
    info.log_path = inquirer.fuzzy(
        message="Which?",
        choices=reversed(info.log_paths)
    ).execute()
    process_log(info)

def quit(info):
    exit()

def make(info, target=None):
    if target is None: 
        target = inquirer.text(message="Target?").execute()
    print_and_run(f"{info.full_snakemake} {target} {info.snakemake_args}")
    update_logs(info)

def select_make(info):
    print(f"{info.snakemake} -Fn")
    info.output_files = re.findall(
        r"    output: (.+)", 
        subprocess.run(
            info.snakemake.split() + ["-Fn"], stdout=subprocess.PIPE
        ).stdout.decode('utf-8')
    )
    targets = inquirer.fuzzy(
        message="Target?", 
        multiselect=True,
        choices=info.output_files
    ).execute()
    make(info, " ".join(targets))


def summarize_value(value):
    if isinstance(value, list): return len(value)
    value = str(value)
    if not value: return value
    return value.splitlines()[0][0:80]

def print_state(info):
    df = pd.DataFrame([
        dict(
            key=key, 
            value=summarize_value(value)
        )
        for key, value in vars(info).items()
    ])
    with pd.option_context('display.max_colwidth', None):
        print(df.sort_values("key"))

def inspect_logs(info):
    log_paths = inquirer.fuzzy(
        message="Which?",
        multiselect=True,
        choices=[
            Choice(row.PATH, name=f"{row.PATH}: {row.STATE}")
            for row in info.jobs_df.itertuples()
        ], 
    ).execute()
    for log_path in log_paths:
        print(f"====== {log_path} ======")
        os.system(f"{info.cmd} {shlex.quote(str(log_path))}")


def interact(): 
    parser = argparse.ArgumentParser(prog='Snakemake interact')
    parser.add_argument("--cmd", default="cat")
    parser.add_argument("--path", default=".")
    parser.add_argument("--hostname", default=socket.gethostname())
    parser.add_argument("--snakemake", default=None)
    parser.add_argument("--screen", default=None)
    info, snakemake_args = parser.parse_known_args()
    info.path = pathlib.Path(info.path)
    info.on_triton = "triton" in info.hostname
    if info.snakemake is None:
        info.snakemake = "poetry run snakemake"
    if info.screen is None:
        info.screen = "screen -dmS 0" if info.on_triton else ""
    info.full_snakemake = f"{info.screen} {info.snakemake}"
    info.snakemake_args = " ".join(map(shlex.quote, snakemake_args)).strip()
    if info.snakemake_args == "auto":
        if info.on_triton:
            info.snakemake_args = "--keep-going --keep-incomplete --slurm -j1024 --default-resources runtime=10 mem_mb=1000 cpus_per_task=1"
        else:
            info.snakemake_args = "-c"

    print_state(info)
    update_logs(info)
    while True:
        print("[Log path]:", info.log_path)
        choices = [update_logs, inspect_logs, make, select_make, inspect_log, select_log, process_log, print_state, quit]
        if info.jobs_df is None: 
            choices.remove(inspect_logs)
        what = inquirer.fuzzy(
            message="What?",
            choices=choices
        ).execute()
        what(info)