import pathlib, subprocess, re, io
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
import pandas as pd
import os


def mtime(path): return pathlib.Path(path).stat().st_mtime

def jobid(path): return int(re.findall(r"([0-9]+)\.log", path)[0])

def path_to_dict(path): return dict(PATH=path, JOBID=jobid(path), STATE="COMPLETED?") | {
    key: "" for key in ["PARTITION", "NAME", "TIME", "START_TIME", "_7"]
}

def update():
    base = pathlib.Path(".snakemake")
    path = sorted((base / "log").glob("*.log"), key=mtime)[-1] 
    return path
    print("Log file:", path)
    with open(path, "r") as fd:
        content = fd.read()
    jobs = re.findall(r"Job .+ has been submitted with SLURM jobid .+ \(log: (.+)\).", content)
    df = pd.DataFrame(map(path_to_dict, jobs) if jobs else [dict(PATH="", STATE="Not ready yet...")])
    for error_path in re.findall(r"log: (.+) \(check log file\(s\) for error details\)", content):
        # print(df.loc[df.path == error_path, "error"])
        df.loc[df.PATH == error_path, "STATE"] = "ERRORED"
    print("# errors:", sum(df.STATE == "ERRORED"))
    progress = re.findall(r"[0-9]+ of [0-9]+ steps \(.+\) done", content)
    if progress: print("progress:", progress[-1])
    # 
    slurm_queue = subprocess.run(['slurm', 'q'], stdout=subprocess.PIPE).stdout.decode('utf-8')
    slurm_df = pd.read_csv(io.StringIO(slurm_queue), sep="\s+")
    # 
    for row in slurm_df.itertuples():
        for key, value in row._asdict().items():
            if key == "Index": continue
            df.loc[df.JOBID == row.JOBID, key] = value
    # 
    print(df[df.STATE != "COMPLETED?"])
    return df

def interact(): 
    cmd = "cat"
    df = update()
    return print(df)
    while True:
        what = inquirer.fuzzy(
            message="What?",
            choices=["print overview", f"change command (current: {cmd})"] + [
                Choice([row.PATH], name=f"{row.PATH}: {row.STATE}")
                for row in df[df.STATE != "COMPLETED?"].itertuples()
            ] + ["print all", "exit"],
        ).execute()
        if what == "print overview": 
            df = update()
            continue
        if what == f"change command (current: {cmd})": 
            cmd = inquirer.text(message=f"Enter command (current: {cmd}):").execute()
            continue
        if what == "exit": break
        if what == "print all":
            what = df.loc[df.STATE != "COMPLETED?", "PATH"]
        df = update()
        for what in what:
            print(f"\n\n======= {what} =======\n\n")
            os.system(f"{cmd} {what}")
            print(f"\n\n======= {what} =======\n\n")