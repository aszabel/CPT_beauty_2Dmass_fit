#include "omp.h"
#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "TH2D.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TChain.h"
#include "TString.h"
#include "TF2.h"
#include "TMath.h"
#include <string>
#include <fstream>
#include "TCanvas.h"
#include "Math/ProbFuncMathCore.h"

typedef std::basic_string<char> string;
typedef std::basic_ifstream<char> ifstream;
typedef std::basic_ofstream<char> ofstream;

using namespace cpt_b0_analysis;
   double minx = 1800.;
   double maxx = 1940.;
   double miny = 2700.;
   double maxy = 8300.;
   const int nvar_md = 7;
   const int nvar_mb = 7;
   const int ncontr = 6;
   const int nbins = 40;

   bool uncorr = false;
   TH2D *hist;

void Draw(const double *res, double nevents, TString name[], int sign);

int main(int argc, char *argv[]){
	if(argc !=2){
	       	std::cout<< " Please state: 'muplus' or 'muminus'." << std::endl;
		return 1;
	}
	if(strcmp(argv[1], "muplus")!=0 && strcmp(argv[1], "mumunus")!=0){
	       	std::cout<< " Please state: 'muplus' or 'muminus'." << std::endl;
		return 1;
	}
	int sign;
	if (strcmp(argv[1], "muplus")==0) sign=1; 
	if (strcmp(argv[1], "muminus")==0) sign=0; 
	
        string minName = "Minuit2";
        string algoName = "";
        ROOT::Math::Minimizer* min =
        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

        // set tolerance , etc...
        min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
        min->SetMaxIterations(10000);  // for GSL
        //min->SetTolerance(2.50e5);
        min->SetTolerance(50.);
        min->SetPrintLevel(2);

        //min->SetStrategy(2);
        //min->SetPrecision(0.00001);




        // create funciton wrapper for minmizer
        // a IMultiGenFunction type
        TChain ch("BlindedTree");
   //ch.Add("/mnt/home/share/lhcb/CPT_beauty/data2016/selected/selected_data2016MagDown.root");
   ch.Add("/home/szabelskia/LHCb/data2016/tree_missPT_D_M_MagDown25102023_nomassDmuCut/selected_data2016MagDown.root");
     double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
   bool charge;


   std::vector<std::pair<double, double>> vect_2D;
   hist = new TH2D("hist", "", nbins, minx, maxx, nbins, miny, maxy);  
   ch.SetBranchAddress("B_M", &B_M);
   ch.SetBranchAddress("missPT", &missPT);
   ch.SetBranchAddress("D_M", &D_M);
   ch.SetBranchAddress("mu_PT", &mu_PT);
   ch.SetBranchAddress("mu_P", &mu_P);
   ch.SetBranchAddress("mu_eta", &mu_eta);
   ch.SetBranchAddress("K_PT", &K_PT);
   ch.SetBranchAddress("truecharge", &charge);


   for (int i=0; i<ch.GetEntries(); ++i){
   //for (int i=0; i<5e5; ++i){
      ch.GetEntry(i);
      if (mu_PT<500 || mu_P<5000 || mu_eta<2 || mu_eta>4.5) continue;
      double B_MMcorr = B_M +2.0*missPT;
      if (B_MMcorr<miny || B_MMcorr>maxy) continue;
      if (D_M<minx || D_M>maxx) continue;
      if (int(charge) == sign){
              vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
              hist->Fill(D_M, B_MMcorr);
      }
   }
        //double  nevents = double(vect_2D.size());
        md_fit md_shape(minx, maxx, ncontr);
	mb_2misspt_fit mb_shape(miny, maxy, ncontr);

	double xx[ncontr][nvar_md];
   	double dxx[ncontr][nvar_md];
   	TString name[ncontr] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm", "sidebands", "Bu2D0Dsm"};
   	for (int ifile=0; ifile<ncontr; ifile++){
		if (ifile==4) continue;
		string filename((TString::Format("D_M_results/res_%s_%d.txt", name[ifile].Data(), sign)).Data());
		std::ifstream infile(filename,  std::ios::binary);
           	int k=0;
           	while(infile>>xx[ifile][k]>> dxx[ifile][k]){
			xx[4][k]=1.0;
			dxx[4][k]=0.001;
			
                	k++;
           	}
   	}
	xx[4][0] = 1.88518e-08;
	dxx[4][0] = 2.97034e-07;
	xx[4][1] = 1.31834e-07;
	dxx[4][1] = 7.13695e-11;

	double xx_mcorr[ncontr][nvar_mb];
	double dxx_mcorr[ncontr][nvar_mb];
	for (int ifile=0; ifile<ncontr; ifile++){
                ifstream infile_mb(TString::Format("B_M_results/res_%s_%d.txt", name[ifile].Data(), sign));
                int k=0;
                while(infile_mb>>xx_mcorr[ifile][k]>> dxx_mcorr[ifile][k]){
                        k++;
                }
        }


        auto fchi2 = [&md_shape, &mb_shape, vect_2D, xx, dxx, xx_mcorr, dxx_mcorr](const double *par)->double{
               	const double *pa = &par[ncontr*(nvar_md+nvar_mb)];
		double frac [ncontr];
		frac[0] = 1-abs(pa[0]);
		frac[1] = abs(pa[0])*(1.0-abs(pa[1]));
    		frac[2] = abs(pa[0])*abs(pa[1])*(1.0-abs(pa[2]));
    		frac[3] = abs(pa[0])*abs(pa[1])*abs(pa[2])*(1.0-abs(pa[3]));
   		frac[4] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4]));
    		frac[5] = abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])*abs(pa[4]);
		double param[ncontr*(nvar_md+nvar_mb)];
		for (int i=0; i<ncontr; i++){
			for (int ivar=0; ivar<nvar_md; ivar++){
				param[i*nvar_md+ivar] = par[i*nvar_md+ivar];
				if (ivar==1 && i!=2 && i!=4)param[i*nvar_md+ivar] = par[1];
			}
		}
	       for (int i=0; i<ncontr; i++){
                        for (int ivar=0; ivar<nvar_mb; ivar++){
				int ipar = ncontr*nvar_md+i*nvar_mb+ivar;
                                param[ipar] = par[ipar];
                        }
                }


		for (int i=0; i<ncontr; i++){
			if (i==4){
				md_shape.int_gaus[i] = 1.0;
				md_shape.int_DCB[i] = 1.0;
				continue;
			}
			double sigma = param[i*nvar_md];
			double mean = param[i*nvar_md+1];
			md_shape.int_gaus.push_back(ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean)); 
			sigma = param[i*nvar_md+2];
			double n = param[i*nvar_md+5];
			double alpha = param[i*nvar_md+4];
			double alpha_h = param[i*nvar_md+6];

			md_shape.int_DCB[i] = TMath::Abs(-ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2.*mean-maxx, alpha_h, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha_h, n, sigma, mean));
		

	                double a1 = par[i*nvar_md+0];
 	                double a2 = par[i*nvar_md+1];


			md_shape.int_Cheb[i] = (1.0-a2)*maxx + 0.5*a1*maxx*maxx+2./3.*a2*maxx*maxx*maxx- (1.0-a2)*minx-0.5*a1*minx*minx-2./3.*a2*minx*minx*minx;

			double s = abs(param[ncontr*nvar_md+i*nvar_mb+1]);
			mean = param[ncontr*nvar_md+i*nvar_mb];
			double minnB = mean-s;
                        if (minx>mean-s) minnB = minx;
                        double maxxB = mean+s;
                        if (maxx<mean+s) maxxB = maxx;

			mb_shape.int_cos[i] = 1.0/(2.0)*(1.0+(maxxB-mean)/s+TMath::Sin((maxxB-mean)/s*TMath::Pi())/TMath::Pi())-1.0/(2.0)*(1.0+(minnB-mean)/s+TMath::Sin((minnB-mean)/s*TMath::Pi())/TMath::Pi());	


			double mean1 = param[ncontr*nvar_md+i*nvar_mb+3];
			double sigma1 = abs(param[ncontr*nvar_md+i*nvar_mb+4]);

			mb_shape.int_gaus1[i] = ROOT::Math::normal_cdf(maxy, sigma1, mean1)-ROOT::Math::normal_cdf(miny, sigma1, mean1);


			double mean2 = param[ncontr*nvar_md+i*nvar_mb+5];
			double sigma2 = abs(param[ncontr*nvar_md+i*nvar_mb+6]);

			mb_shape.int_gaus2[i] = ROOT::Math::normal_cdf(maxy, sigma2, mean2)-ROOT::Math::normal_cdf(miny, sigma2, mean2);
		}

 		double chi2 = 0.0;
	        #pragma omp parallel for reduction (+:chi2)
                for (auto &v: vect_2D){
			double mdass = std::get<0>(v);
			double mcorr = std::get<1>(v);
			double likelihood = 0.0;
			for (int i=0; i<ncontr; i++){
				double md_like, mb_like;
				md_like = md_shape.func_full(&mdass, &param[i*nvar_md], i);
				mb_like = mb_shape.func_full(&mcorr, &param[ncontr*nvar_md+i*nvar_mb], i);
				likelihood+= md_like*mb_like*frac[i];
			}
                        if(likelihood>0.0) chi2 -= 2.0*log(likelihood);
		}
		for (int i=0; i<ncontr; i++){
			if(i==4) continue;
			for (int ivar=0; ivar<nvar_md; ivar++){
				if(dxx[i][ivar]!=0)chi2+= (xx[i][ivar]-param[i*nvar_md+ivar])*(xx[i][ivar]-param[i*nvar_md+ivar])/(dxx[i][ivar]*dxx[i][ivar]);   // use the results of MC fits
                	}
		}
		for (int i=0; i<ncontr; i++){
                        for (int ivar=0; ivar<nvar_mb; ivar++){
                        	if(dxx_mcorr[i][ivar]!=0)chi2+= (xx_mcorr[i][ivar]-param[ncontr*nvar_md+i*nvar_mb+ivar])*(xx_mcorr[i][ivar]-param[ncontr*nvar_md+i*nvar_mb+ivar])/(dxx_mcorr[i][ivar]*dxx_mcorr[i][ivar]);   // use the results of MC fits
                        }
                }
		
                return chi2;
        };


	TString varname_md[nvar_md] = {"sigma", "mean", "sigmaCB", "f12", "alpha", "n", "alpha_h"};

        for (int i=0; i<ncontr; i++){
		for (int ivar=0; ivar<nvar_md; ivar++){
			min->SetVariable(i*nvar_md+ivar, (TString::Format("%s_%s", name[i].Data(), varname_md[ivar].Data())).Data(), xx[i][ivar], dxx[i][ivar]+1.0e-11);
			//min->FixVariable(i*nvar_md+ivar);
			min->SetVariableLowerLimit(i*nvar_md+ivar, 1e-13);
			if(ivar==3)min->SetVariableLimits(i*nvar_md+ivar, -1.0, 1.0);

			if (ivar==5 || ivar==6) min->FixVariable(i*nvar_md+ivar);
			if (i==4 && ivar!=0 && ivar!=1) min->FixVariable(i*nvar_md+ivar); // combinatorial
			//if (i==4) min->FixVariable(i*nvar_md+ivar); // combinatorial
			//if(ivar!=1&&ivar!=0)min->FixVariable(i*nvar_md+ivar);
			if(i!=0 && ivar==1) min->FixVariable(i*nvar_md+ivar);
		}
	}


	
	
	
	TString varname_mb[nvar_mb] = {"mean_rc", "sigma_rc", "f12_gaus", "mean1_gaus", "sigma1_gaus", "mean2_gaus", "sigma2_gaus"};

        for (int i=0; i<ncontr; i++){
                for (int ivar=0; ivar<nvar_mb; ivar++){
                        min->SetVariable(ncontr*nvar_md+i*nvar_mb+ivar, (TString::Format("%s_%s", name[i].Data(), varname_mb[ivar].Data())).Data(), xx_mcorr[i][ivar], dxx_mcorr[i][ivar]+1.0e-11);
			//min->FixVariable(ncontr*nvar_md+i*nvar_mb+ivar);
                        min->SetVariableLowerLimit(ncontr*nvar_md+i*nvar_mb+ivar, 1e-13);
                        if(ivar==2)min->SetVariableLimits(ncontr*nvar_md+i*nvar_mb+ivar, -1.0, 1.0);

                }
        }



      	//      for (auto i: vect_fix_mb) min->FixVariable(ncontr*nvar_md+i);
        double frac_init[ncontr-1] = {0.468636, 0.736568, 0.973902, 0.949009, -0.00476566};
	double frac_init_plus[ncontr-1] = {0.578309, 0.816153, 0.519927, 0.627151, 0.0151929}; 
        //double frac_init[ncontr-1] = {0.657945, 0.901491, 0.968305, 0.527389, 0.0034851};
        for (int i=0; i<ncontr-1; i++){

                if (sign==0) min->SetVariable((nvar_md+nvar_mb)*ncontr+i, (TString::Format("par_frac%d", i)).Data(), frac_init[i], 0.01);
		else min->SetVariable((nvar_md+nvar_mb)*ncontr+i, (TString::Format("par_frac%d", i)).Data(), frac_init_plus[i], 0.01);

                min->SetVariableLimits((nvar_md+nvar_mb)*ncontr+i, -1.0, 1.0);
                //min->FixVariable((nvar_md+nvar_mb)*ncontr+i);
        }
	ROOT::Math::Functor f(fchi2, (nvar_md+nvar_mb)*ncontr+ncontr-1);
	min->SetFunction(f);
	double CL_normal = ROOT::Math::normal_cdf(1) -  ROOT::Math::normal_cdf(-1);

        min->SetErrorDef(TMath::ChisquareQuantile(CL_normal, min->NFree()));

	double x_fix[ncontr*(nvar_md+nvar_mb)+ncontr-1];
	double dx_fix[ncontr*(nvar_md+nvar_mb)+ncontr-1];
	ifstream input_fix(TString::Format("results_fix_%d.txt", sign));
	int k=0;
	while (input_fix >> x_fix[k] >> dx_fix[k]) k++; 
	input_fix.close();
	uncorr = true;
	if (k==ncontr*(nvar_md+nvar_mb)+ncontr-1) std::cout << "TUUUUUUU \n\n";

	/*
 	for(int ivar=0; ivar<ncontr*(nvar_md+nvar_mb)+ncontr-1; ivar++){ 
		std::cout << x_fix[ivar] << "  " << ivar << std::endl;
		min->SetVariableValue(ivar, x_fix[ivar]);
	}*/
        min->Minimize();
	bool over[ncontr*(nvar_md+nvar_mb)];
	for (int i=0; i<ncontr*(nvar_md+nvar_mb); i++){
		over[i]=false;
		for (int j=0; j<ncontr-1; j++){
			double cor = min->Correlation(i, ncontr*(nvar_md+nvar_mb)+j);
			if (abs(cor) < 0.02) {
				std::cout << min->VariableName(i) << "  " << cor << std::endl;
				
			}else over[i]= true;
		}
	}
	std::vector<int> vect_fix = {5, 8, 12, 15, 19, 22, 26, 29, 30, 31, 32, 33, 34, 36, 40, 57, 62};
	for (int i=0; i<ncontr*(nvar_md+nvar_mb)+ncontr-1; i++) min->SetVariableValue(i, x_fix[i]);
	uncorr = true;	
	for (auto i: vect_fix) if (uncorr) min->FixVariable(i);

	//min->Minimize();

	ofstream results(TString::Format("results_%d.txt", sign));
	for (int i=0; i<(nvar_md+nvar_mb)*ncontr+ncontr-1; i++){
		results<< min->X()[i]<< "  " << min->Errors()[i] << std::endl;
	}
	results.close();

    return 0;	   
}
