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

varname = "varname_md"
directory = "MC_directory_MD"
fig_prefix = "Dmass"
shapes = "DMshapes"

scan_path = (Path("scans") / ("scan_sidebands_" + config_path.stem + ("_binned" if config["binned"] else "_unbinned"))).absolute()
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
    for sign in ["muminus", "muplus"]:
        for i, fit in enumerate(config["contrName"]):
            for j, var in enumerate(config[varname][i]):
                limits = config["scanLimitsVect"][f"{fit}_{var}"]
                v = random.uniform(limits[0], limits[1])
                config["init_values"][i][j] = v

        config["sign"] = sign
        out_path = scan_path / str(n) / sign
        out_path.mkdir(parents=True, exist_ok=True)
        out_name = out_path / (config_path.stem + ".json")
        with out_name.open("w") as f:
            json.dump(config, f, indent=2)

        os.chdir(out_path)
        with open(".rootrc", "w") as f:
            f.write(f"Unix.*.Root.MacroPath:    .:{scripts_path.parent}/macros:$(ROOTSYS)/macros")
        run_all = scripts_path.parent / "fit1D_mass"
        print(f"Running fit {n} ...")
        p = Popen([str(run_all), str(out_name.name)], stdout=PIPE, stderr=PIPE)
        outs, errs = p.communicate()
        with open("out.log", "wb") as f:
            f.write(outs)
        with open("err.log", "wb") as f:
            f.write(errs)

        result_paths = list(Path(".").glob("res_sidebands_*.txt"))
        if not result_paths:
            print("Missing results directory - FIT FAILED")
            continue
        with result_paths[0].open() as f:
            data = f.readlines()
            print(f"Status={data[0].strip().split()[0]}, Chi2={data[0].strip().split()[1]}")

results = {}
for run in scan_path.iterdir():
    if run.name == "results":
        continue
    for sign_path in run.iterdir():
        sign = sign_path.name
        result_paths = list(sign_path.glob("res_sidebands_*.txt"))
        if not result_paths:
            continue
        fname = result_paths[0]
        with fname.open() as f:
            data = f.readlines()
            status = int(data[0].strip().split()[0])
            if status <= 1:
                if sign in results:
                    v = float(data[0].strip().split()[1])
                    if v < results[sign]["best_value"]:
                        results[sign]["best_run"] = run.stem
                        results[sign]["best_value"] = v
                    if math.fabs(v - 1.0) < math.fabs(results[sign]["good_value"] - 1.0):
                        results[sign]["good_run"] = run.stem
                        results[sign]["good_value"] = v
                else:
                    results[sign] = {
                        "best_run": run.stem,
                        "best_value": float(data[0].strip().split()[1]),
                        "good_run": run.stem,
                        "good_value": float(data[0].strip().split()[1])
                    }

for sign in results:
    print(f"Best run for {sign}: ", results[sign]["best_run"])
    output_path = scan_path / "results" / sign
    best_run_path = scan_path / results[sign]["best_run"] / sign
    shutil.copytree(best_run_path, output_path)

for run in scan_path.iterdir():
    if run.name == "results":
        continue
    shutil.rmtree(run)

pprint(results)
