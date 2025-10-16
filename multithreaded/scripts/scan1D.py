#!/bin/env python3
import sys
import os
import json
import random
import math
import shutil
import argparse

from pathlib import Path
from subprocess import Popen, PIPE
from pprint import pprint

scripts_path = Path( __file__ ).absolute().parent

parser = argparse.ArgumentParser(
    description='Random search of 1D mass fit input parameters'
)
parser.add_argument('-c', '--config', type=str, required=True, help="Config file")
parser.add_argument('-n', '--tries', type=int, default=100, help="Number of random search tries")
parser.add_argument('-i', '--input', type=str, help="Path with fit result txt files to use as random search range")

args = parser.parse_args()

config_path = Path(args.config)
N = args.tries
with config_path.open() as f:
    config = json.load(f)

if config["category"] == "1D_BM":
    varname = "varname_mb"
    script = "runBMall.sh"
    directory = "MC_directory_MB"
    fig_prefix = "Bmass"
    shapes = "BMshapes"
elif config["category"] == "1D_DM":
    varname = "varname_md"
    script = "runDMall.sh"
    directory = "MC_directory_MD"
    fig_prefix = "Dmass"
    shapes = "DMshapes"
else:
    raise(Exception("Unknown fit category: " + config["category"]))

scan_path = (Path("scans") / ("scan_" + config_path.stem + ("_binned" if config["binned"] else "_unbinned"))).absolute()
scan_path.mkdir(parents=True, exist_ok=True)

if args.input:
    input_path = Path(args.input)
    for results_path in input_path.iterdir():
        tags = results_path.stem.split("_")
        fit = tags[1]
        pdf = tags[2]
        idx = config["contrName"].index(fit)
        assert pdf == config[shapes][idx], f"Non compatible pdf parameters: {pdf} vs {config[shapes][idx]}"
        with results_path.open() as results:
            data = results.readlines()
            for j, var in enumerate(config[varname][idx]):
                mean, sigma = (float(v) for v in data[j+2].strip().split())
                config["scanLimitsVect"][f"{fit}_{var}"] = [
                    mean - sigma, mean + sigma
                ]

for n in range(N):
    for i, fit in enumerate(config["contrName"]):
        for j, var in enumerate(config[varname][i]):
            limits = config["scanLimitsVect"][f"{fit}_{var}"]
            v = random.uniform(limits[0], limits[1])
            config["init_values"][i][j] = v

    out_path = scan_path / str(n)
    out_path.mkdir(parents=True, exist_ok=True)
    out_name = out_path / (config_path.stem + ".json")
    with out_name.open("w") as f:
        json.dump(config, f, indent=2)

    os.chdir(out_path)
    with open(".rootrc", "w") as f:
        f.write(f"Unix.*.Root.MacroPath:    .:{scripts_path.parent}/macros:$(ROOTSYS)/macros")
    run_all = scripts_path / script
    print(f"Running fit {n} ...")
    p = Popen([str(run_all), str(out_name.name)], stdout=PIPE, stderr=PIPE)
    outs, errs = p.communicate()
    with open("out.log", "wb") as f:
        f.write(outs)
    with open("err.log", "wb") as f:
        f.write(errs)

    result_path = Path("results") / (config[directory] + ("_binned" if config["binned"] else "_unbinned"))
    if not result_path.exists():
        print("Missing results directory - FIT FAILED")
        continue
    for fname in result_path.iterdir():
        with fname.open() as f:
            data = f.readlines()
            print(f"{fname.stem} : Status={data[0].strip().split()[0]}, Chi2={data[0].strip().split()[1]}")

results = {}
for run in scan_path.iterdir():
    if run.name == "results":
        continue
    result_path = run / "results" / (config[directory] + ("_binned" if config["binned"] else "_unbinned"))
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
                    if math.fabs(v - 1.0) < math.fabs(results[fname.stem]["good_value"] - 1.0):
                        results[fname.stem]["good_run"] = run.stem
                        results[fname.stem]["good_value"] = v
                else:
                    results[fname.stem] = {
                        "best_run": run.stem,
                        "best_value": float(data[0].strip().split()[1]),
                        "good_run": run.stem,
                        "good_value": float(data[0].strip().split()[1])
                    }

output_path = scan_path / "results"
output_path.mkdir(parents=True, exist_ok=True)
for fit, data in results.items():
    fig = fit.replace("res", fig_prefix)
    best_fit_path = output_path / fit / "best"
    best_fit_path.mkdir(parents=True, exist_ok=True)
    best_run_path = scan_path / data["best_run"]
    best_result_path = best_run_path / "results" / (config[directory] + ("_binned" if config["binned"] else "_unbinned"))
    best_figures_path = best_run_path / "results" / (config[directory] + ("_binned" if config["binned"] else "_unbinned") + "_figures")
    for f in best_run_path.iterdir():
        if not f.is_dir():
            shutil.copy(f, best_fit_path / f.name)
    shutil.copy(best_result_path / (fit + ".txt"), best_fit_path / (fit + ".txt"))
    shutil.copy(best_figures_path / (fig + ".pdf"), best_fit_path / (fig + ".pdf"))

    good_fit_path = output_path / fit / "good"
    good_fit_path.mkdir(parents=True, exist_ok=True)
    good_run_path = scan_path / data["good_run"]
    good_result_path = good_run_path / "results" / (config[directory] + ("_binned" if config["binned"] else "_unbinned"))
    good_figures_path = good_run_path / "results" / (config[directory] + ("_binned" if config["binned"] else "_unbinned") + "_figures")
    for f in good_run_path.iterdir():
        if not f.is_dir():
            shutil.copy(f, good_fit_path / f.name)
    shutil.copy(good_result_path / (fit + ".txt"), good_fit_path / (fit + ".txt"))
    shutil.copy(good_figures_path / (fig + ".pdf"), good_fit_path / (fig + ".pdf"))

for run in scan_path.iterdir():
    if run.name == "results":
        continue
    shutil.rmtree(run)

pprint(results)
