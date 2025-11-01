# CPT_beauty_2Dmass_fit

The fit startegy:
- Obtain the shapes from MC and m_D sidebands using 1D fits
- Obtain the combinatorial background fraction estimate from 1D fit to m_D
- Obtain a conservative estimate of contribution fractions from 2D fit with shapes **fixed** to MC values and using the combinatorial background fraction as a starting point
- Run the the fit starting from previous ones with shape parametres constrained to MC values

Each step can be repeated with randomly smeared input values ...

## 1D fits

At first we perform 1D fits to MC and data sidebands to extract expected shapes of various contributions. Those shapes will be used as constraints in the final 2D fit.

1. File `../common/src/BasicShapes.cpp` contains all relevant PDF components.
2. Files `common/src/M_B_2missPT_fit.cpp` and `common/src/D_M_fit_shape.cpp` contain the shapes of the `m_Bcorr` and `m_D` repectively.
3. Files `macros/B_M_fit_Every.C` and `macros/D_M_fit_Every.C` are used to perform fits to all contributions
4. Files `configs/config_1D_*.json` contain fit parameters as well as paths to the MC on shared CIS space.
5. The fit is executed with `./scripts/runBMall.sh configs/config_1D_BM.json` or `./scripts/runDMall.sh configs/config_1D_DM.json`
Running the above produces new `results` subdirectories named according to setup stored in the config file. For each fit two directories are created `B_M_<config>_results_<binned>` (`D_M_<config>_results_<binned>`) and `B_M_<config>_results_<binned>_figures` (`D_M_<config>_results_<binned>_figures`)
The values in results are used in the final fits as starting points and gaussian constraints.
The text files with fit results contain: the fit status and the Chi2/NLL in the firts line, parameter values with associated error each on separate line.
6. A random search of the fit input parameters can be performed. The results will be stored in the `scans` directory.

Either provide the parameter scan ranges in the config file via:
```
"scanLimitsVect":{
  "signal_sigma": [100, 3000],
  "signal_mean": [1000, 10000],
  "signal_alpha": [0.1, 5.0],
  "signal_n": [1, 1],
  "signal_xi": [1000, 10000]
}
```

To execute a random search of input parameters with 1000 tries use:
```
./sbatch-scan1D.sh configs/config_1D_BM_nominal.json 1000
```

Or use results of a previous fit to set the ranges according to parameter uncertainties.

To execute a random search of input parameters with 1000 tries for "signal" only with limits on random search set to 2.0 sigma from a previously executed fit run:
```
./sbatch-scan1D.sh configs/config_1D_BM_nominal.json 1000 limits/nominal_BM signal 2.0
```

## Sidebands fit

The first estimate of the combinatorial background fraction is extracted from a fit to sidebands of the m_D.

1. The fit is executed using the `fit1D_mass` executable: `./fit1D_mass configs/config_sidebands.json`
2. The last line in the `res_sidebands_0.txt` or `res_sidebands_1.txt` files contains the sidebands fraction and its error.
3. A random search of the fit input parameters can be performed. The results will be stored in the `scans` directory.

To execute a random search of input parameters with 1000 tries use:
```
./sbatch-sidebands.sh configs/config_sidebands.json 1000
```

## 2D fits

1. **??? DEPRECATED ???** With the code you can run 1D fits in the mcorr and m_D dimentions: `root -l 'run_Dfit_DM.C(sign)'` or `root -l 'run_Dfit_BM.C(sign)'`
2. **DEPRECATED ???** The 2D fit can be run with `root -l 'run_Dfit_2D.C(sign)'`, in the `macros/Dfit_2D.C` the 2D-pdfs are defined as multiplications of the two 1D-pdfs.

1. The fit is executed using the `fit2D_mass` binary
2. There are several fit options available:
   - "frac" - a fit with floating fraction and all parameters fixed to the values from 1D mass fits
   - "all" - a fit with floating fraction and all parameters constrained to the values from 1D mass fits
   - other options fix either m_B or m_D
3. The final results are stored in results_0.txt and results_1.txt
4. Draw_details.C produces 2D and 1D plots using data (for the moment the path is set to my PC location) and from the results stored in one of the above files.
5. The Draw_details.C macro can by run using `./runDraw.sh configs/config_2D_nominal.json`
6. The fit is repeated until a valid one is found, subsequent fits have input variables randomly smeared
7. The fit config file allows to set a number of tries - this results in additional fit repetitions with the best selected as a final results
8. It is possible to schedule several fits in parallel with different random seeds using the `./sbatch-scan2D.sh` script

To schedule the fractions fit 100 times to slurm use:
```
./sbatch-scan2D.sh configs/config_2D_nominal_fractions.json 100
```

To schedule the final fit 100 times to slurm use:
```
./sbatch-scan2D.sh configs/config_2D_nominal_fractions.json 100
```

## sWeights

TODO

## LifeTime fit

TODO

# Bibliography

## Likelihood and Simultanous fits

* https://indico.belle2.org/event/3456/contributions/18541/attachments/10210/15687/fitting-belle2-academy.pdf
* https://arxiv.org/pdf/physics/0401045v1
* https://stat.ufl.edu/wp-content/uploads/sites/120/voneshTalk.pdf
* https://www.researchgate.net/publication/12410534_Parameter_correlations_while_curve_fitting
* https://root-forum.cern.ch/t/errors-and-contours-of-fit-parameters/41373/5

# TODO

3. Computational Approach: Principal Component Analysis (PCA)

If you are performing numerical optimization for parameter estimation (like MLE or method of moments), you can use a computational method to decorrelate the final parameter space:

    Initial Fit: Perform an initial fit to obtain the estimated parameters θ=(γ,δ,ξ,λ).

    Estimate Covariance: Calculate the covariance matrix Σ for the estimated parameters. This is often the inverse of the observed or expected Fisher Information Matrix.

    PCA/Eigen-decomposition: Apply Principal Component Analysis (PCA) or perform an eigen-decomposition on the covariance matrix Σ:
    Σ=VΛVT

    where V is the matrix of eigenvectors (principal components) and Λ is the diagonal matrix of eigenvalues (variances).

    Transform Parameters: Define the new, uncorrelated parameters η as a linear transformation of the original parameters θ:
    η=VTθ

The new parameters η will have a diagonal covariance matrix Λ, meaning they are uncorrelated. You can then try to re-optimize or analyze the results in this decorrelated η space, which can sometimes aid convergence and interpretation.

ROOT code to get eigenvectors and eigenvalues
https://root-forum.cern.ch/t/eigenvalues-of-a-nearly-singular-covariance-matrix/18694/5
