# CPT_beauty_2Dmass_fit

1. Files `common/src/M_B_2missPT_fit.cpp` and `common/src/D_M_fit_shape.cpp` contain the shapes of the mcorr (raised_cosine+2gaus) and m_D (gaus + doublesidedcrystalball) repectively.
2. files `macros/B_M_fit_Every.C` and `macros/D_M_fit_Every.C` contain paths to the MC on shared CIS space for MagDown. They are run with `root -l 'run_B_M_Every.C(choice, sign)'` ( `'Burun_B_M_Every.C(choice, sign)'`) where choice is between B0->Dmunu(signal), Bu->Dmunu, Bs->Dsmunu, B0->DsD, combinatorial(sidebands), Bu->DsD. 'sign' is the sign of the muon.
Running the above produces two directories B_M_results (D_M_results) and B_M_figures(D_M_figures). The valuse in results are used in the final fits as starting points and gaussian constraints. To run all the MC shape fits you can run
`./runBMall.sh` and `./runDMall.sh`
3. With the code you can run 1D fits in the mcorr and m_D dimentions: `root -l 'run_Dfit_DM.C(sign)'` or `root -l 'run_Dfit_BM.C(sign)'`
4. The 2D fit can be run with `root -l 'run_Dfit_2D.C(sign)'`, in the `macros/Dfit_2D.C` the 2D-pdfs are defined as multiplications of the two 1D-pdfs.
5. The final results are stored in results_0.txt and results_1.txt
6. Draw_details.C produces 2D and 1D plots using data (for the moment the path is set to my PC location) and from the results stored in one of the above files.
