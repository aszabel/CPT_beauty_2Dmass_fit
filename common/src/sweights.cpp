#include "sweights.h"
#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "config.h"

#include "TROOT.h"
#include "TChain.h"
#include "TMatrixD.h"
#include "TDecompLU.h" 

#include <vector>
#include <iostream>

#include "TFile.h"

namespace cpt_b0_analysis
{
	sWeights::sWeights(const std::string& _path_data, const std::string& _outfilename, const std::string& _outtreename): path_data(_path_data), outfilename(_outfilename), outtreename(_outtreename){};
        
	void  sWeights::get_sWeigths(const double res[], bool BBbar){	
	std::vector<double> vec_Tau, vec_sWeights;
	std::string plus_minus;
	if (BBbar) plus_minus = "plus";
	else plus_minus = "minus";
	TFile *outfile = new TFile(Form("%s_%s.root",outfilename.c_str(), plus_minus.c_str()), "recreate");
       	TTree *outtree = new TTree(outtreename.c_str(), outtreename.c_str());
	outtree -> Branch("vec_Tau", &vec_Tau);	
	outtree -> Branch("vec_sWeights", &vec_sWeights);	
		int ncontr = Config::ncontr;
		int nvar_md = Config::nvar_md;
		int nvar_mb = Config::nvar_mb;
		double frac[ncontr];
		for (int i=0; i<ncontr; i++){
			frac[i] = std::abs(res[ncontr*(nvar_md+nvar_mb)+i]);
			std::cout << frac[i] << " frac " << i << std::endl;
		}

		TChain ch(Config::chainName.c_str());

		double Tau, D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
  		bool charge;

		ch.Add(path_data.c_str());
		std::vector<std::pair<double, double>> vect_2D={};
   		ch.SetBranchAddress("Tau", &Tau);
   		ch.SetBranchAddress("B_M", &B_M);
   		ch.SetBranchAddress("missPT", &missPT);
   		ch.SetBranchAddress("D_M", &D_M);
   		ch.SetBranchAddress("mu_PT", &mu_PT);
   		ch.SetBranchAddress("mu_P", &mu_P);
   		ch.SetBranchAddress("mu_eta", &mu_eta);
   		ch.SetBranchAddress("K_PT", &K_PT);
   		ch.SetBranchAddress("truecharge", &charge);
   		//ch.SetBranchAddress("blinded_charge", &bcharge);

		std::vector<double> vect_time = {};
		std::cout << " tuuuuuu\n";

		int Nevents = Config::nentries;
		if (Nevents==-1)
			Nevents = ch.GetEntries();
		for (int i=0; i<Nevents;i++){
      			ch.GetEntry(i);
			if (mu_PT < Config::muPTmin || mu_P < Config::muPmin || mu_eta < Config::eta_min || mu_eta > Config::eta_max)
                                continue;
      			double B_MMcorr = B_M +2.0*missPT;
      			if (B_MMcorr < Config::minBMcorr || B_MMcorr > Config::maxBMcorr)
                        continue;

                	if (D_M < Config::minDM || D_M > Config::maxDM) continue;
			if (Tau<Config::tMin || Tau>Config::tMax) continue; 
      			if (charge==BBbar){
				vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
				vect_time.push_back(Tau);
			}
			
      		}
		std::cout << vect_2D.size() << " siiiizieeeee \n";

		TMatrixD matrixV(ncontr, ncontr);
		const auto& B_PDFs = Config::getVectorPDFs("Bmass");
		const auto& D_PDFs = Config::getVectorPDFs("Dmass");
		double nevents = (double) vect_2D.size();
		for (int i=0; i<ncontr; i++){
			D_PDFs[i]->CalcIntegral(&res[i * Config::nvar_md], Config::minDM, Config::maxDM);
                        B_PDFs[i]->CalcIntegral(&res[Config::ncontr * Config::nvar_md + i * Config::nvar_mb], Config::minBMcorr, Config::maxBMcorr);
		}
		for (int n=0; n<ncontr; n++){
			for(int j=0; j<ncontr; j++){
				matrixV (n, j)=0.0;
				int emax = (int)vect_2D.size();
	        		for (int e=0; e<emax; e++){
					double mdass = std::get<0>(vect_2D[e]);
                        		double mcorr = std::get<1>(vect_2D[e]);

					double sum_k = 0.0;
					double md_like, mb_like, md_like2, mb_like2;
		    			for (int i=0; i<ncontr; i++){
						md_like = D_PDFs[i]->EvalPDF(&mdass, &res[i * Config::nvar_md]);
                	                	mb_like = B_PDFs[i]->EvalPDF(&mcorr, &res[Config::ncontr * Config::nvar_md + i * Config::nvar_mb]);
				//std::cout << md_like << "   " << mb_like << "  " << nevents << std::endl;
						sum_k += md_like*mb_like*frac[i]*nevents;	
					}
					//std::cout << sum_k << " ssssssummmmmmm\n";
					md_like = D_PDFs[n]->EvalPDF(&mdass, &res[n * Config::nvar_md]);
					mb_like = B_PDFs[n]->EvalPDF(&mcorr, &res[Config::ncontr * Config::nvar_md + n * Config::nvar_mb]);
					md_like2 = D_PDFs[j]->EvalPDF(&mdass, &res[j * Config::nvar_md]);
					mb_like2 = B_PDFs[j]->EvalPDF(&mcorr, &res[Config::ncontr * Config::nvar_md + j * Config::nvar_mb]);
					matrixV (n, j) += md_like*mb_like*md_like2*mb_like2/(sum_k*sum_k);  
				}
			}
	      	}	
   		TDecompLU lu(matrixV);
		bool invert=false;
		TMatrixD inv_matrixV = lu.Invert(invert);
		std::vector<std::pair<double, double>> sWeights={};
		for (int e=0; e<nevents; e++){
                        double mdass = std::get<0>(vect_2D[e]);
                        double mcorr = std::get<1>(vect_2D[e]);
                        double Tau = vect_time[e];

			double sum_k = 0.0;
                        double md_like, mb_like;
                        for (int i=0; i<ncontr; i++){
                        	md_like = D_PDFs[i]->EvalPDF(&mdass, &res[i * Config::nvar_md]);
                                mb_like = B_PDFs[i]->EvalPDF(&mcorr, &res[Config::ncontr * Config::nvar_md + i * Config::nvar_mb]);
                                sum_k += md_like*mb_like*frac[i]*nevents;
                        }

			double sum_n = 0.0;
			for (int i=0; i<ncontr; i++){
				md_like = D_PDFs[i]->EvalPDF(&mdass, &res[i * Config::nvar_md]);
                                mb_like = B_PDFs[i]->EvalPDF(&mcorr, &res[Config::ncontr * Config::nvar_md + i * Config::nvar_mb]);
				sum_n+=inv_matrixV(0, i)*md_like*mb_like;
			}
			vec_Tau.push_back(Tau);
			vec_sWeights.push_back(sum_n/sum_k);
			sWeights.push_back(std::make_pair(Tau, sum_n/sum_k));
		//	std::cout<< sum_n << "  " << sum_k << std::endl;
		}
		outtree->Fill();
		outtree->Write();
		outfile->Close();
		std::cout << vec_sWeights.size() << " siiiizieeeee " << vec_Tau.size() << "\n";
	}	

}
