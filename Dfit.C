#include "Dmass_fit.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


using namespace cpt_b0_analysis;

const int nbins = 40;
   //const int ncontr = 6;   

void Dfit(){

    std::string minName = "Minuit2";
 std::string algoName = "";
 ROOT::Math::Minimizer* min =
      ROOT::Math::Factory::CreateMinimizer(minName, algoName);

// set tolerance , etc...
   min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
   min->SetMaxIterations(10000);  // for GSL
   min->SetTolerance(5.0e3);
   min->SetPrintLevel(2);

   double sigma_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);
   int nfree = 7;

   min->SetErrorDef(ROOT::Math::chisquared_quantile(sigma_normal2D, nfree));
   //min->SetStrategy(2);
   //min->SetPrecision(0.00001);

   // create funciton wrapper for minmizer
   // a IMultiGenFunction type

   TChain ch("BlindedTree");
   //ch.Add("/home/szabel/LHCb/data2016/asl_selection_tree/tree_missPT_D_MM/a_sl_selected_ntuple*.root");
   //ch.Add("/home/szabelskia/LHCb/Semileptonic_analysis_project/Semileptonic_analysis_project/strip28r1MagUp2016/newSelection/simpleDmass2/full_selection_newTau_allCuts07_D_M_PT_rap_2016MagUp.root");
   //ch.Add("/home/szabelskia/LHCb/Semileptonic_analysis_project/Semileptonic_analysis_project/strip28r1MagUp2016/newSelection/simpleDmass2/full_selection_newTau_allCuts07_D_M_PT_rap_2016MagDown.root");
   ch.Add("/mnt/home/aszabelski/data2016/strip28r1_Selection/tree_missPT_D_M_MagDown25102023_nomassDmuCut/selected_data2016MagDown.root");

   //TFile fileh("hDp_mass_mu_250bins_range1800_1940_22016magUp.root", "recreate");
   //TFile fileh("hDp_mass_mu_250bins_range1800_1940_2016magUp_tight_eta.root", "recreate");
   double minx = 1800.;
   double maxx = 1940.;
   double miny = 2700.;
   double maxy = 8300.;
   std::vector<std::pair<double, double>> vect_Dmass;
   
   double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
   bool charge;

   TH2D *hcorr[ncontr];
   hcorr[0] = new TH2D("hcorr0", "", 2, minx, maxx, nbins, miny, maxy);
   ch.SetBranchAddress("B_M", &B_M);
   ch.SetBranchAddress("missPT", &missPT);
   ch.SetBranchAddress("D_M", &D_M);
   ch.SetBranchAddress("mu_PT", &mu_PT);
   ch.SetBranchAddress("mu_P", &mu_P);
   ch.SetBranchAddress("mu_eta", &mu_eta);
   ch.SetBranchAddress("K_PT", &K_PT);
   ch.SetBranchAddress("truecharge", &charge);


   bool sign = false;
   for (int i=0; i<ch.GetEntries(); ++i){
      ch.GetEntry(i);
      if (mu_PT<500 || mu_P<5000 || mu_eta<2 || mu_eta>4.5) continue;
      double B_MMcorr = sqrt(B_M * B_M) +2.0*TMath::Abs(missPT);
      if (B_MMcorr<miny || B_MMcorr>maxy) continue;
      if (D_M<minx || D_M>maxx) continue;
      if (charge == sign){
	      vect_Dmass.push_back(std::make_pair(D_M, B_MMcorr));
	      hcorr[0]->Fill(D_M, B_MMcorr);
      }
   }
   for (int i=1; i<ncontr; i++){
	   hcorr[i] = (TH2D*)hcorr[0]->Clone(Form("hcorr%d", i));
   }
   
   cout << vect_Dmass.size() << " entriesssss \n";

   double xx[ncontr][nprm];
   double dxx[ncontr][nprm];
   TString name[ncontr] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm" , "sidebands", "Bu2D0Dsm"};
   for (int ifile=0; ifile<ncontr; ifile++){
           ifstream infile(Form("res_%s_0.txt", name[ifile].Data()));
           int k=0;
           while(infile>>xx[ifile][k]>> dxx[ifile][k]){
                cout << xx[ifile][k] << "  " << dxx[ifile][k] << endl;
                k++;
           }
   }
   const int nsigBu=6;
   TString name_sigBu[nsigBu] = {"sigma_sigBu", "mean_sigBu", "f12_sigBu", "sigmaCB_sigBu", "alpha_sigBu", "n_sigBu"};
   double data_sigBu[nsigBu];
   double err_sigBu[nsigBu];
   ifstream input_sigBu("res_D_M_sig_Bu_0.txt");
           int k=0;
           while(input_sigBu>>data_sigBu[k]>> err_sigBu[k]){
                k++;
           }
   const int nBs=4;
   TString name_Bs[nBs] = {  "sigmaCB_Bs", "mean_Bs", "alpha_Bs", "n_CB_Bs"};

   double data_Bs[nBs];
   double err_Bs[nBs];
   ifstream input_Bs("res_D_M_Bs_0.txt");
           k=0;
           while(input_Bs>>data_Bs[k]>> err_Bs[k]){
                k++;
           }
   const int nB0D=7;
   TString name_B0D[nB0D] = {"sigma_B0D", "mean_B0D", "sigmaCB_B0D", "f12_B0D", "alpha_B0D", "n_B0D", "alphah_B0D"};

   double data_B0D[nB0D];
   double err_B0D[nB0D];
   ifstream input_B0D("res_D_M_B0D_0.txt");
           k=0;
           while(input_B0D>>data_B0D[k]>> err_B0D[k]){
                k++;
           }
   const int nBuD=7;
   TString name_BuD[nBuD] = {"sigma_BuD", "mean_BuD", "sigmaCB_BuD", "f12_BuD", "alpha_BuD", "n_BuD", "alphah_BuD"};

   double data_BuD[nBuD];
   double err_BuD[nBuD];
   ifstream input_BuD("res_D_M_BuD_0.txt");
           k=0;
           while(input_BuD>>data_BuD[k]>> err_BuD[k]){
                k++;
           }
   double yy[nsigBu+ nBs+ nB0D+nBuD];
   double dyy[nsigBu+ nBs+ nB0D+nBuD];

   for (int isigBu=0; isigBu<nsigBu; isigBu++){
	   yy[isigBu] = data_sigBu[isigBu];
	   dyy[isigBu] = err_sigBu[isigBu];
   }
   for (int iBs=0; iBs<nBs; iBs++){
	   yy[nsigBu+iBs] = data_Bs[iBs];
	   dyy[nsigBu+iBs] = err_Bs[iBs];
   }
   for (int iB0D=0; iB0D<nB0D; iB0D++){
	   yy[nsigBu+nBs+iB0D] = data_B0D[iB0D];
	   dyy[nsigBu+nBs+iB0D] = err_B0D[iB0D];
   }
   for (int iBuD=0; iBuD<nBuD; iBuD++){
	   yy[nsigBu+nBs+nB0D+iBuD] = data_BuD[iBuD];
	   dyy[nsigBu+nBs+nB0D+iBuD] = err_BuD[iBuD];
   }
 

   TFile *file[ncontr];
   TTree *t[ncontr];
   const int npar = 8;
   double alpha = 2.40868;
   double nCB = 0.5;
   double sigma = 9.65705;
   double mean = 1869.49;
   double exp_slope = 0.00477;
   double s12 = 5.96338;
   double f12 = 0.475104;
   double fcomb = 0.250746;
 
   const double par[npar] = {sigma, mean, exp_slope, s12, f12, fcomb, alpha, nCB};


   Dmass Dfit(vect_Dmass, hcorr, minx, maxx, miny, maxy, sign);
   Dfit.SetExternal(xx, dxx, yy, dyy);

   ROOT::Math::Functor f(Dfit, 87);
  


   min->SetFunction(f);



   for (int i=0; i<nsigBu; i++){
	min->SetVariable(i, name_sigBu[i].Data(), data_sigBu[i], err_sigBu[i]+1.0e-6);
	//min->FixVariable(i);
	if (i==5) min->FixVariable(i);
	if (i==4) min->SetVariableLimits(i, 0.1, 100.);	
   }	

   for (int i=0; i<nBs; i++){
	min->SetVariable(nsigBu+i, name_Bs[i].Data(), data_Bs[i], err_Bs[i]+1.0e-6);
	//min->FixVariable(nsigBu+i);
	if (i==3) min->FixVariable(nsigBu+i);
	if (i==2) min->SetVariableLimits(nsigBu+i, 0.1, 100.);	
   }	

   for (int i=0; i<nB0D; i++){
	min->SetVariable(nsigBu+nBs+i, name_B0D[i].Data(), data_B0D[i], err_B0D[i]+1.0e-6);
	//min->FixVariable(nsigBu+nBs+i);
	if (i==5 || i==1 || i==6) min->FixVariable(nsigBu+nBs+i);
	if (i==4) min->SetVariableLimits(nsigBu+nBs+i, 0.1, 100.);	
   }	

   for (int i=0; i<nBuD; i++){
	min->SetVariable(nsigBu+nBs+nB0D+i, name_BuD[i].Data(), data_BuD[i], err_BuD[i]+1.0e-6);
	//min->FixVariable(nsigBu+nBs+nB0D+i);
	if (i==5 || i==1 || i==6) min->FixVariable(nsigBu+nBs+nB0D+i);
	if (i==4) min->SetVariableLimits(nsigBu+nBs+nB0D+i, 0.1, 100.);	
   }	






   min->SetVariable(24, "exp_slope", 0.00477, 0.0001);
   min->FixVariable(24);//

   min->SetVariable(25, "a1", -0.00102419, 0.000001);
   //min->FixVariable(25);
   
   min->SetVariable(26, "a2", 1.31649e-07, 1.e-9);
   //min->FixVariable(26);

   min->SetVariable(27, "nevents", 150000., 100.);
   min->FixVariable(27);


   TString var[nprm] = {"mean", "sigma", "norm", "gnorm1", "gmean1", "gsigma1", "gnorm2", "gmean2", "gsigma2"};

   double prev_min[nprm*ncontr] = {
  5887.77, 1047.44, 1, 1, 5515.27, 900.9, 0.207176, 4974.8, 424.313,
5331.12, 1948.47, 1, 1, 3850.68, 571.252, 0.802861, 5208.66, 293.803,
5683.08, 1939.56, 1, 1, 4857.33, 704.641, 0.798811, 5930.25, 355.965,
4650.74, 2196.15, 1, 1, 4574.9, 840.197, 0.323682, 7217.14, 1278.25,
5508.03, 1545.6, 1, 1, 6637.3, 975.379, 0.725038, 4197.69, 685.763,
4661.07, 2206.7, 1, 1, 4568.66, 844.516, 0.185074, 5945.48, 1243.97};
   for (int i=0; i<ncontr; i++){
           for (int j=0; j<nprm; j++){
                   min->SetVariable(28+i*nprm+j, Form("%s_%s", name[i].Data(), var[j].Data()), prev_min[i*nprm+j], dxx[i][j]+1.e-6);
                   if(dxx[i][j]==0) {
                            min->FixVariable(28+i*nprm+j);
                   }
	   }
   }
        const std::string name_frac[ncontr-1] ={"fBu", "fCB", "fDCB", "fcomb", "fD"};

   double frac_par_fitted[ncontr-1]={0.629472, 0.847221, 0.5915, 0.598838, 0.300802};
   for(int ic=0; ic<ncontr-1; ic++){
           min->SetVariable(82+ic,name_frac[ic],frac_par_fitted[ic], 0.001);
   }

   

   min->SetVariableLimits(82, 0.0, 1.0);

   min->SetVariableLimits(83, 0.0, 1.0);


   min->SetVariableLimits(84, 0.0, 1.0);

   min->SetVariableLimits(85, 0.0, 1.0);

   min->SetVariableLimits(86, 0.0, 1.0);

   double pars[87], errors[87];

   ifstream input("results_500k_fixed_frac.txt");
   k=0;
   while (input>> pars[k] >> errors[k])k++;	

   for(int ipar = 0; ipar<87; ipar++) min->SetVariableValue(ipar, pars[ipar]);
  
	

 
   //for (int i=0; i<82; i++) if(i!=25||i!=26) min->FixVariable(i); 

   //min->FixVariable(25);	
   //min->FixVariable(26);	
   min->Minimize();

   double CovMtx[nvar*nvar];
   min->GetCovMatrix(CovMtx);

   Dfit.Draw(min->X(), "Dmass", CovMtx);

   const double *res = min->X();
   const double *err = min->Errors();

   ofstream output("results.txt");
   for (int i=0; i<nvar; i++){ 
   	output << res[i] << "  " << err[i] << endl;
   }
   output.close();
}


