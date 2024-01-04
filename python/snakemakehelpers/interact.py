import pathlib, subprocess, re, io
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
import pandas as pd
import os, shlex
import argparse
from argparse import Namespace
import socket
import shutil
import random, string

def random_string(length, letters=string.ascii_letters):
    return ''.join(random.choice(letters) for i in range(length))

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
    info.jobs_df = None
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
    info.jobs_df = None
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
    print(info.jobs_df.loc[
        info.jobs_df.STATE.isin(info.states) if info.states else info.jobs_df.STATE != "",
        info.columns
    ])
    print(info.jobs_df.groupby(["STATE"])["STATE"].count())

    
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
    if info.screen:
        print("Started screen session in background, refresh logs manually in a bit...")
    update_logs(info)

def select_make(info):
    print(f"{info.snakemake} -Fn")
    process = subprocess.run(
        info.snakemake.split() + ["-Fn"], 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout = process.stdout.decode('utf-8')
    if process.returncode != 0:
        return print(f"Snakemake errored:\n\n{stdout}")
    info.output_files = re.findall(
        r"    output: (.+)", 
        stdout
    )
    if not info.output_files:
        return print("No possible output files found!")

    targets = inquirer.fuzzy(
        message="Target?", 
        multiselect=True,
        choices=info.output_files
    ).execute()
    make(info, " ".join(targets))


def summarize_value(value):
    if isinstance(value, list): 
        return len(value)
    value = str(value)
    if not value: 
        return value
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

def modify_columns(info):
    info.columns = inquirer.fuzzy(
        message="Which?",
        multiselect=True,
        choices=info.jobs_df.columns
    ).execute()

def modify_states(info):
    info.states = inquirer.fuzzy(
        message="Which?",
        multiselect=True,
        choices=sorted(info.jobs_df.STATE.unique())
    ).execute()

def interact(): 
    parser = argparse.ArgumentParser(prog='Snakemake interact')
    parser.add_argument("--cmd", default="cat")
    parser.add_argument("--path", default=".")
    parser.add_argument("--hostname", default=socket.gethostname())
    parser.add_argument("--slurm", default=None)
    parser.add_argument("--snakemake", default=None)
    parser.add_argument("--screen", default=None)
    parser.add_argument("--jobs", default=256)
    parser.add_argument("--runtime", default=10)
    parser.add_argument("--mem", default=1000)
    parser.add_argument("--cpu", default=1)
    parser.add_argument("--snakemake-args", default="auto")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-c", "--local", action="store_true")
    parser.add_argument("targets", nargs="*")
    info = parser.parse_args()
    info.path = pathlib.Path(info.path)
    info.on_triton = "triton" in info.hostname
    if info.slurm is None:
        if shutil.which("slurm") and not info.local and not info.dry_run: info.slurm = True
        else: info.slurm = False
    if info.snakemake is None:
        if shutil.which("snakemake"): info.snakemake = "snakemake"
        elif shutil.which("poetry"): info.snakemake = "poetry run snakemake"
        else: raise ValueError("Could not determine snakemake executable")
    if info.screen is None:
        if info.slurm and shutil.which("screen"): info.screen = "screen -dmS " + random_string(4)
        else: info.screen = "" 
    info.full_snakemake = f"{info.screen} {info.snakemake}"
    # info.snakemake_args = " ".join(map(shlex.quote, snakemake_args)).strip()
    if info.snakemake_args in ["auto", ""]:
        if info.slurm:
            info.snakemake_args = f"--keep-going --keep-incomplete --slurm -j{info.jobs} --default-resources runtime={info.runtime} mem_mb={info.mem} cpus_per_task={info.cpu}"
        else:
            info.snakemake_args = "-c"
            if info.dry_run: 
                info.snakemake_args = "-n"
    info.columns = ["PATH", "STATE", "TIME"]
    info.states = []

    if info.targets:
        make(info, " ".join(map(shlex.quote, info.targets)))
    print_state(info)
    update_logs(info)
    while True:
        print("[Log path]:", info.log_path)
        if getattr(info, "progress", None) is not None: print(info.progress)
        choices = [
            update_logs, inspect_logs, select_make, 
            modify_columns, modify_states, 
            make, inspect_log, select_log, print_state, 
            quit
        ]
        if info.jobs_df is None: 
            choices.remove(inspect_logs)
            choices.remove(modify_columns)
            choices.remove(modify_states)
        descriptions = dict(
            update_logs=f"Process latest snakemake log file",
            inspect_logs=f"Inspect (multiple) job log files ({info.cmd} ...)",
            select_make=f"Select target(s) to make ({info.snakemake} -Fn)",
            make=f"Make target(s) ({info.full_snakemake} ... {info.snakemake_args})",
            inspect_log=f"Inspect selected snakemake log file ({info.log_path})",
            modify_columns=f"Modify columns of jobs that get printed ({info.columns})",
            modify_states=f"Modify states of jobs that get printed ({info.states})",
            select_log=f"Select a different snakemake log file",
            process_log=f"Process selected snakemake log file ({info.log_path})",
            print_state=f"Print internal state of this program",
            quit="quit"
        )
        choices = [
            Choice(func, name=descriptions.get(func.__name__)) for func in choices
        ]
        what = inquirer.fuzzy(
            message="What?",
            choices=choices
        ).execute()
        what(info)