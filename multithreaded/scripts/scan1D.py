#!/bin/env python3
import sys
import os
import json
import random
import math
import shutil

from pathlib import Path
from subprocess import Popen, PIPE
from pprint import pprint

scripts_path = Path( __file__ ).absolute().parent

config_path = Path(sys.argv[1])
N = int(sys.argv[2])
with config_path.open() as f:
    config = json.load(f)

if config["category"] == "1D_BM":
    varname = "varname_mb"
    script = "runBMall.sh"
    directory = "MC_directory_MB"
    fig_prefix = "Bmass"
elif config["category"] == "1D_DM":
    varname = "varname_md"
    script = "runDMall.sh"
    directory = "MC_directory_MD"
    fig_prefix = "Dmass"
else:
    raise(Exception("Unknown fit category: " + config["category"])



scan_path = (Path("scans") / ("scan_" + config_path.stem + "_binned" if config["binned"] else "_unbinned")).absolute()
scan_path.mkdir(parents=True, exist_ok=True)

for n in range(N):
    for i, fit in enumerate(config["Fits"]):
        for j, var in enumerate(config[varname]):
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

    result_path = Path("results") / (config[directory] + "_binned" if config["binned"] else "_unbinned")
    for fname in result_path.iterdir():
        with fname.open() as f:
            data = f.readlines()
            print(f"{fname.stem} : Status={data[0].strip()}, Chi2={data[1].strip()}")

results = {}
for run in scan_path.iterdir():
    if run.name == "results":
        continue
    result_path = run / "results" / (config[directory] + "_binned" if config["binned"] else "_unbinned")
    for fname in result_path.iterdir():
        with fname.open() as f:
            data = f.readlines()
            status = int(data[0].strip())
            if status <= 1:
                if fname.stem in results:
                    v = float(data[1].strip())
                    if v < results[fname.stem]["best_value"]:
                        results[fname.stem]["best_run"] = run.stem
                        results[fname.stem]["best_value"] = v
                    if math.fabs(v - 1.0) < math.fabs(results[fname.stem]["good_value"] - 1.0):
                        results[fname.stem]["good_run"] = run.stem
                        results[fname.stem]["good_value"] = v
                else:
                    results[fname.stem] = {
                        "best_run": run.stem,
                        "best_value": float(data[1].strip()),
                        "good_run": run.stem,
                        "good_value": float(data[1].strip())
                    }

output_path = scan_path / "results"
output_path.mkdir(parents=True, exist_ok=True)
for fit, data in results.items():
    fig = fit.replace("res", fig_prefix)
    best_fit_path = output_path / fit / "best"
    best_fit_path.mkdir(parents=True, exist_ok=True)
    best_run_path = scan_path / data["best_run"]
    best_result_path = best_run_path / "results" / (config[directory] + "_binned" if config["binned"] else "_unbinned")
    best_figures_path = best_run_path / "results" / (config[directory] + ("_binned" if config["binned"] else "_unbinned") + "_figures")
    for f in best_run_path.iterdir():
        if not f.is_dir():
            shutil.copy(f, best_fit_path / f.name)
    shutil.copy(best_result_path / (fit + ".txt"), best_fit_path / (fit + ".txt"))
    shutil.copy(best_figures_path / (fig + ".pdf"), best_fit_path / (fig + ".pdf"))

    good_fit_path = output_path / fit / "good"
    good_fit_path.mkdir(parents=True, exist_ok=True)
    good_run_path = scan_path / data["good_run"]
    good_result_path = good_run_path / "results" / (config[directory] + "_binned" if config["binned"] else "_unbinned")
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
