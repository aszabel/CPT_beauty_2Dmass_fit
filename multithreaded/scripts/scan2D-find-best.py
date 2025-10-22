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
args = parser.parse_args()

config_path = Path(args.config)
with config_path.open() as f:
    config = json.load(f)

scan_path = (Path("scans2D") / ("scan_" + config_path.stem + ("_binned" if config["binned"] else "_unbinned"))).absolute()
print("Analysing: ", scan_path)

results = {}
chi2 = []
for run in scan_path.iterdir():
    if not str(run.name).startswith("fit2D"):
        continue
    result_path = run / ("best_results" + ("_binned" if config["binned"] else "_unbinned")) / "fit2D_best"
    if not result_path.exists():
        continue
    for fname in result_path.iterdir():
        with fname.open() as f:
            data = f.readlines()
            status = int(data[0].strip().split()[0])
            if status <= 1:
                if fname.stem in results:
                    v = float(data[0].strip().split()[1])
                    if v < results[fname.stem]["best_value"]:
                        results[fname.stem]["best_run"] = run.stem
                        results[fname.stem]["best_value"] = v
                else:
                    results[fname.stem] = {
                        "best_run": run.stem,
                        "best_value": float(data[0].strip().split()[1]),
                    }

output_path = scan_path / ("best_results" + ("_binned" if config["binned"] else "_unbinned"))
figure_path = scan_path / ("best_results" + ("_binned" if config["binned"] else "_unbinned") + "_figures")
(output_path / "fit2D_best").mkdir(parents=True, exist_ok=True)
figure_path.mkdir(parents=True, exist_ok=True)
#with (output_path / "chi2.json").open("w") as f:
#    json.dump(chi2, f)

with open(output_path.parent / ".rootrc", "w") as f:
    f.write(f"Unix.*.Root.MacroPath:    .:{scripts_path.parent}/macros:$(ROOTSYS)/macros")

for fit, data in results.items():
    sign = "muminus" if fit.split("_")[1] == "0" else "muplus"
    best_run_path = scan_path / data["best_run"]
    best_result_path = best_run_path / ("best_results" + ("_binned" if config["binned"] else "_unbinned")) / "fit2D_best"
    os.chdir(best_run_path)

    best_config = ""
    for f in best_run_path.glob(f"config*{sign}*"):
        if not f.is_dir():
            if str(f).endswith(sign):
                best_config = f
                with f.open() as data_file:
                    data = json.load(data_file)
                mc_path = Path(data["MC_directory_MD"])
                data["MC_directory_MD"] = str(mc_path.resolve())
                mc_path = Path(data["MC_directory_MB"])
                data["MC_directory_MB"] = str(mc_path.resolve())
                with (output_path / f.name).open("w") as out_file:
                    json.dump(data, out_file, indent=2)
            else:
                shutil.copy(f, output_path / f.name)
    shutil.copy(best_result_path / (fit + ".txt"), output_path / "fit2D_best" / (fit + ".txt"))

    (figure_path / sign).mkdir(parents=True, exist_ok=True)
    subprocess.check_call(
        [
            scripts_path / "runDraw.sh",
            output_path / best_config.name,
            output_path.parent,
            figure_path / sign
        ]
    )

pprint(results)