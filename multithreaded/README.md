# CPT_beauty_2Dmass_fit

## 1D fits

At first we perform 1D fits to MC and data sidebands to extract expected shapes of various contributions. Those shapes will be used as constraints in the final 2D fit.

1. File `../common/src/BasicShapes.cpp` contains all relevant PDF components.
2. Files `common/src/M_B_2missPT_fit.cpp` and `common/src/D_M_fit_shape.cpp` contain the shapes of the `m_Bcorr` and `m_D` repectively.
3. Files `macros/B_M_fit_Every.C` and `macros/D_M_fit_Every.C` are used to perform fits to all contributions
4. Files `configs/config_1D_*.json` contain fit parameters as well as paths to the MC on shared CIS space.
5. The fit is executed with `./runBMall.sh configs/config_1D_BM.json` or `./runDMall.sh configs/config_1D_DM.json`
Running the above produces new `results` subdirectories named according to setup stored in the config file. For each fit two directories are created `B_M_<config>_results_<binned>` (`D_M_<config>_results_<binned>`) and `B_M_<config>_results_<binned>_figures` (`D_M_<config>_results_<binned>_figures`)
The values in results are used in the final fits as starting points and gaussian constraints.

To execute a random search of input parameters with 1000 tries use:
```
sbatch ./scripts/scan1D.slurm configs/config_1D_BM.json 1000
```

## 2D fits

1. With the code you can run 1D fits in the mcorr and m_D dimentions: `root -l 'run_Dfit_DM.C(sign)'` or `root -l 'run_Dfit_BM.C(sign)'`
2. The 2D fit can be run with `root -l 'run_Dfit_2D.C(sign)'`, in the `macros/Dfit_2D.C` the 2D-pdfs are defined as multiplications of the two 1D-pdfs.
3. The final results are stored in results_0.txt and results_1.txt
4. Draw_details.C produces 2D and 1D plots using data (for the moment the path is set to my PC location) and from the results stored in one of the above files.
