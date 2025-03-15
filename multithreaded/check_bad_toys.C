#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
using namespace std;



void check_fits(int &best, int itoy, string results_dir, string results_num){
	double min_chi2 = 1.0e15;
	double chi2;
	int status;
	best = -1;
 	int nfits = 100;
	if (results_dir == "results_binned100kMU_2") nfits=30;
	for (int i=1; i<=nfits; i++){

		ostringstream os;
	 	os << "toy_res/toy_" << itoy << "/" << results_dir.c_str()<<"/fit2D_" << i << "/results_" << results_num.c_str() << ".txt";
		string str1 = os.str();
		ifstream resin(str1);
		resin>>status >> chi2;
		resin.close();
		if ((status==0)&& chi2<min_chi2){
			min_chi2 = chi2;
			best = i;
		}
	}
	if (best==-1 || min_chi2>24.8e6) cout << itoy << " ";
}
void check_bad_toys(){

	gStyle->SetOptFit(1);
	const int nfrac = 6;
	const int nres = 2;
		double res[(7+7)*(nfrac)+2*nfrac];
		double dres[(7+7)*(nfrac)+2*nfrac];
	TCanvas *c[nfrac][nres];
        for (int ifrac=0; ifrac<nfrac; ifrac++){
		for (int ires=0; ires<nres; ires++){
                	c[ifrac][ires] = new TCanvas(Form("c%d_%d", ifrac, ires), "", 500, 500);
		}
	}  
	TH1D *pull[nfrac][nres];
	double gen_frac[nfrac] = {0.58, 0.06, 0.027, 0.06, 0.27, 0.0}; 
	const int ntoys = 500;
	for (int ifrac=0; ifrac<nfrac; ifrac++){
		for (int ires=0; ires<nres; ires++){
			ostringstream os;
	 		os << "pull" << ifrac << "_" << ires;
			string str1 = os.str();
			pull[ifrac][ires] = new TH1D( str1.c_str(), "", 200, gen_frac[ifrac]-0.03, gen_frac[ifrac]+0.03);
		}
	}

	string results_num[nres] = {"0_0", "0_3"};
	//string results_num[nres] = {"0_0", "0_2", "0_3", "0_3"};
	//string results_dir[nres] = {"results_binned100kMU_2", "results_binned100kMU_2", "results_unbinned100kMU_2"}; 
	string results_dir[nres] = {"results_binned100kMU_2",  "results_unbinned100kMU_2"}; 

  //      TFile file("results_tuple.root", "recreate");
	int status[nres] ={-1, -1};//, -1, -1};
	double chi2[nres]= {0.0, 0.0};// , 0.0, 0.0};
	double resfrac[nres][nfrac];
	int toynum[nres]={-1, -1};//, -1, -1};
	int best[nres]={-1, -1};//, -1, -1};

	TTree *tree[nres];
	for (int i=0; i<nres; i++){
		ostringstream os;
	 	os << "restree" << i;
		string str1 = os.str();
		tree[i] = new TTree(str1.c_str(), str1.c_str());
		tree[i]->Branch("toynum", &toynum[i], "toynum/I");
		tree[i]->Branch("best", &best[i],  "best/I");
		tree[i]->Branch("status", &status[i], "status/I");
		tree[i]->Branch("chi2", &chi2[i], "chi2/D");
		for (int j=0; j<nfrac; j++){
			ostringstream os;
                	os << "frac" << j;
                	string str1 = os.str();

			tree[i]->Branch(str1.c_str(), &resfrac[i][j], (str1+"/D").c_str());
		}
	}
		


//#pragma omp parallel for	
	for (int itoy=1; itoy<=ntoys; itoy++){
		//cout << "============= \n" << itoy << endl << endl;
		
		for (int ires=0; ires<nres; ires++){
			int bestie=0;
			check_fits(bestie, itoy, results_dir[ires], results_num[ires]);
			ostringstream os;
	 		os << "toy_res/toy_" << itoy << "/" << results_dir[ires].c_str()<<"/fit2D_" << bestie << "/results_" << results_num[ires].c_str() << ".txt";
			string str1 = os.str();
			best[ires] = bestie;
			//ifstream input(Form("toy_res/toy_%d/%s/fit2D_%d/results_%s.txt", itoy, results_dir[ires].c_str(), best, results_num[ires].c_str()));
			ifstream input(str1.c_str());
			
			double x, dx;
                	input>>status[ires] >> chi2[ires];
			int i=0;
        		while (input>>x>>dx){
                 		res[i] = abs(x);
				dres[i] = dx;  
				i++;
        		}
        		input.close();

			int n0 = (nfrac) * (7 + 7);
		
			//cout << "======================================\n";
			for (int ifrac=0; ifrac<nfrac; ifrac++){
					resfrac[ires][ifrac] = res[n0+ifrac];
					pull[ifrac][ires]->Fill(res[n0+ifrac]);
					//cout << res[n0+ifrac] << " frac" << ifrac << " " << dres[n0+ifrac] << endl << endl;
			
			}
			//cout << res[n0] << " frac0\n";
			toynum[ires]=itoy;
			tree[ires]->Fill();
		}
	}
//	tree[0]->Write();
//	tree[1]->Write();
	//tree[2]->Write();
//	file.Close();
	
	for (int ifrac=0; ifrac<nfrac; ifrac++){
		for (int ires=0; ires<nres; ires++){
			c[ifrac][ires]->cd();
			pull[ifrac][ires]->Draw();
			pull[ifrac][ires]->Fit("gaus");
			ostringstream os;
	 		os << "plots/frac" << ifrac << "_" << results_dir[ires].c_str()<<"_" << results_num[ires].c_str() << ".pdf";
			string str1 = os.str();
			//c[ifrac][ires]->SaveAs(TString::Form("plots/frac%d_%s_%s.pdf", ifrac, results_dir[ires].c_str(), results_num[ires].c_str()));
			c[ifrac][ires]->SaveAs(str1.c_str());
		}
	}
}	
