#!/bin/env python3
import sys
import os
import json
import random

from pathlib import Path
from subprocess import Popen, PIPE
from pprint import pprint

script_path = Path( __file__ ).absolute()

config_path = Path(sys.argv[1])
N = int(sys.argv[2])
with config_path.open() as f:
    config = json.load(f)

scan_path = (Path("scans") / ("scan_" + config_path.stem + "_binned" if config["binned"] else "_unbinned")).absolute()
scan_path.mkdir(parents=True, exist_ok=True)

for n in range(N):
    for i, fit in enumerate(config["Fits"]):
        for j, var in enumerate(config["varname_mb"]):
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
        f.write(f"Unix.*.Root.MacroPath:    .:{script_path.parent}/macros:$(ROOTSYS)/macros")
    run_all = script_path.parent / "runBMall.sh"
    print(f"Running fit {n} ...")
    p = Popen([str(run_all), str(out_name.name)], stdout=PIPE, stderr=PIPE)
    outs, errs = p.communicate()
    with open("out.log", "wb") as f:
        f.write(outs)
    with open("err.log", "wb") as f:
        f.write(errs)

    result_path = Path("results") / (config["MC_directory_MB"] + "_binned" if config["binned"] else "_unbinned")
    for fname in result_path.iterdir():
        with fname.open() as f:
            data = f.readlines()
            print(f"{fname.stem} : Status={data[0].strip()}, Chi2={data[1].strip()}")

results = {}
for run in scan_path.iterdir():
    result_path = run / "results" / (config["MC_directory_MB"] + "_binned" if config["binned"] else "_unbinned")
    for fname in result_path.iterdir():
        with fname.open() as f:
            data = f.readlines()
            status = int(data[0].strip())
            if status <= 1:
                if fname.stem in results:
                    v = float(data[1].strip())
                    if v < results[fname.stem]["value"]:
                        results[fname.stem] = {
                            "run": run.stem,
                            "value": v
                        }
                else:
                    results[fname.stem] = {
                        "run": run.stem,
                        "value": float(data[1].strip())
                    }

pprint(results)
