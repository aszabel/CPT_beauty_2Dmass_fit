#!/bin/env python3
import sys
import os
import json
import random
import math
import shutil
import argparse
import subprocess

from pathlib import Path
from pprint import pprint

scripts_path = Path( __file__ ).absolute().parent
logs_path = Path(os.getcwd()) / "logs"

parser = argparse.ArgumentParser(
    description='Random search of 1D mass fit input parameters'
)
parser.add_argument('-c', '--config', type=str, required=True, help="Config file")
parser.add_argument('-n', '--tries', type=int, default=100, help="Number of random search tries")
parser.add_argument('-s', '--start', type=int, default=1, help="Start jobs from this number")

args = parser.parse_args()

script = "scan2D.slurm"

config_path = Path(args.config)
N = args.tries
START = args.start
with config_path.open() as f:
    config = json.load(f)

scan_path = (Path("scans2D") / ("scan_" + config_path.stem + ("_binned" if config["binned"] else "_unbinned"))).absolute()
scan_path.mkdir(parents=True, exist_ok=True)

for n in range(START, N+1):
    config["randSeed"] = n

    out_path = scan_path / f"fit2D_{n}"
    out_path.mkdir(parents=True, exist_ok=True)
    out_log_path = out_path / "logs"
    out_log_path.symlink_to(logs_path)
    out_name = out_path / (config_path.stem + ".json")
    with out_name.open("w") as f:
        json.dump(config, f, indent=2)

    os.chdir(out_path)
    with open(".rootrc", "w") as f:
        f.write(f"Unix.*.Root.MacroPath:    .:{scripts_path.parent}/macros:$(ROOTSYS)/macros")
    print(f"Scheduling fit {n} ...")
    subprocess.check_call(
        [
            "sbatch",
            "-p",
            "INTEL_HASWELL,INTEL_CASCADE,INTEL_SKYLAKE",
            "-J",
            f"scan2D-{config_path.stem}",
            scripts_path / script,
            scripts_path.parent,
            out_name
        ]
    )