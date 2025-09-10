# CPT_beauty_2Dmass_fit

## 1D fits

At first we perform 1D fits to MC and data sidebands to extract expected shapes of various contributions. Those shapes will be used as constraints in the final 2D fit.

1. Files `common/src/M_B_2missPT_fit.cpp` and `common/src/D_M_fit_shape.cpp` contain the shapes of the mcorr (raised_cosine+2gaus) and m_D (gaus + doublesidedcrystalball) repectively.
2. Files `macros/B_M_fit_Every.C` and `macros/D_M_fit_Every.C` are used to perform fits to all contributions
3. Files `configs/config_1D_*.json` contain fit parameters as well as paths to the MC on shared CIS space.
4. The fit is executed with `./runBMall.sh configs/config_1D_BM.json` or `./runDMall.sh configs/config_1D_DM.json`
Running the above produces two directories `B_M_results` (`D_M_results`) and `B_M_figures` (`D_M_figures`).
The values in results are used in the final fits as starting points and gaussian constraints.

## 2D fits

1. With the code you can run 1D fits in the mcorr and m_D dimentions: `root -l 'run_Dfit_DM.C(sign)'` or `root -l 'run_Dfit_BM.C(sign)'`
2. The 2D fit can be run with `root -l 'run_Dfit_2D.C(sign)'`, in the `macros/Dfit_2D.C` the 2D-pdfs are defined as multiplications of the two 1D-pdfs.
3. The final results are stored in results_0.txt and results_1.txt
4. Draw_details.C produces 2D and 1D plots using data (for the moment the path is set to my PC location) and from the results stored in one of the above files.
