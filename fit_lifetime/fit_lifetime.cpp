#include "pdfs_cpt_model.h"
#include "decay_parameters.h"
#include "sweights.h"
#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "ChebyshevPDF.h"


#include <string>
#include <fstream>
#include <iomanip>
#include "omp.h"
#include "json.hpp"
#include "config_time.h"

// ROOT includes
#include "TROOT.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/ProbFuncMathCore.h"
#include "FastSum.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"

#include <functional>
#include <memory>
#include <tuple>

#include "TF1.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"

#include <mutex>
std::mutex my_mutex;

// bool avx = false;
bool avx = true;


using json = nlohmann::json;
using namespace cpt_b0_analysis;
std::function<double(const double*)> wrap_chi2(TH1D *logHist_expectB, TH1D *logHist_expectBbar, TH1D *logAsym, TH1D *logExpAsym, std::function<double(double, const double*, bool)> func_conv, std::vector<double> vec_time, std::vector<double> vec_time_bar,  TH1D *hist_B, TH1D *hist_Bbar, TH1D *logHistB, TH1D *logHistBbar, bool make_asym, std::vector<std::pair<int, int>> replaceIndexVect);

TH1D* rebin_to_log_bins(TH1D *histin, char *name );

int main(int argc, char *argv[]){
	
	
	// The first, fundamental operation to be performed in order to make ROOT
        // thread-aware.
        ROOT::EnableThreadSafety();

       try{
                if (argc != 2)
                {
                        std::cerr << "Usage: fit2D_mass config.json" << std::endl;
                        return 1;
                }
                if (Config::load(argv[1])){
                        std::cerr<< " Bad config file! " << std::endl;
                        return 1;
                }

        }
        catch (const std::exception& ex) {
                std::cerr << "Error: " << ex.what() << std::endl;
                return 1;
        }




	pdfs PDFs;
	auto func_f = [&PDFs] (double *x, double *par) -> double{ return PDFs.func1D_f(x, par);};

	TF1 pdf_f("pdf_f", func_f, Params::tMin, Params::tMax, 4);
	pdf_f.SetParameters(Params::Rez, 0.2, Params::Gamma, Params::eta);
	std::cout << pdf_f.Integral(Params::tMin, Params::tMax) << std::endl;	

	auto func_fbar = [&PDFs] (double *x, double *par) -> double{ return PDFs.func1D_fbar(x, par);};

	TF1 pdf_fbar("pdf_fbar", func_fbar, Params::tMin, Params::tMax, 4);
	pdf_fbar.SetParameters(Params::Rez, Params::Imz, Params::Gamma, Params::eta);
	std::cout << pdf_fbar.Integral(Params::tMin, Params::tMax) << std::endl;	
	
	auto func_acc = [&PDFs] (double *x, double *par) -> double{ return PDFs.time_acceptance(x, par);};

	TF1 pdf_acc("pdf_acc", func_acc, Params::tMin, Params::tMax, 3);
	pdf_acc.SetParameters(Params::t_shift, Params::alpha, Params::beta);
	std::cout << pdf_acc.Integral(Params::tMin, Params::tMax) << std::endl;	


	TFile file_r(Config::r_factor_file.c_str(), "readonly");
	auto tree_r = (TTree*)file_r.Get("BlindedTree");
	bool charge;
	double truetau, tau;
	tree_r->SetBranchAddress("blinded_charge", &charge);
	tree_r->SetBranchAddress("Tau", &tau);
	tree_r->SetBranchAddress("TrueTau", &truetau);

	gSystem->Exec("mkdir -p trees");
	TFile filetree(Form("trees/tree_hists_%d.root", Config::randSeed), "recreate");
	//TFile filedraw(Form("dawings/draw_hists_%d.root", Config::randSeed), "recreate");
	int Nrbins = 100;
	TH1D hist_rB("hist_rB", "", Nrbins, Config::r_factor_min, Config::r_factor_max);
	TH1D hist_rBbar("hist_rBbar", "", Nrbins, Config::r_factor_min, Config::r_factor_max);

	double *logBins = new double [Config::nlogBins+1];

    // Logarithmic binning: Calculate the bin edges
    for (int i = 0; i <= Config::nlogBins; ++i) {
        logBins[i] = TMath::Power(10, (i * (TMath::Log10(Params::tMax) - TMath::Log10(Params::tMin)) / Config::nlogBins) + TMath::Log10(Params::tMin));
    }

    
    	TH1D *logHistB = new TH1D("logHistB", "", Config::nlogBins, logBins);
    	TH1D *logHistBbar = new TH1D("logHistBbar", "", Config::nlogBins, logBins);
    	TH1D *logHist_expectB = new TH1D("logHist_resB", "", Config::nlogBins, logBins);
    	TH1D *logHist_expectBbar = new TH1D("logHist_resBbar", "", Config::nlogBins, logBins);	
	TH1D *hist_B = new TH1D("hist_B", "", Config::Nbins, Params::tMin, Params::tMax);
	TH1D *hist_Bbar = new TH1D("hist_Bbar", "", Config::Nbins, Params::tMin, Params::tMax);

	//fitted polynomial of truetau vs tau
	double Apol = -0.0114834;
        double Bpol = 0.82811; 
        double Cpol = 0.00869728; 

	std::vector<double> vec_time = {};
	std::vector<double> vec_time_bar = {};
	for (int i=0; i<tree_r->GetEntries(); i++){
		tree_r->GetEntry(i);
		if (tau<Params::tMin || tau>Params::tMax) 
			continue;
		if (charge){
		        hist_rB.Fill((truetau)/(Apol*(tau)*(tau)+Bpol*tau+Cpol));	
		}else{
		       	hist_rBbar.Fill((truetau)/(Apol*tau*tau+Bpol*tau+Cpol));
		}
	}

	hist_rB.Scale(1./hist_rB.Integral());
	hist_rBbar.Scale(1./hist_rBbar.Integral());
	hist_rB.Write();
	hist_rBbar.Write();
	std::vector<double> r_bin_center, r_bin_cont, r_bin_center_bar, r_bin_cont_bar;
	for (int r_bin=1; r_bin<=Nrbins; ++r_bin){
     		r_bin_center.push_back(hist_rB.GetBinCenter(r_bin));
     		r_bin_cont.push_back(hist_rB.GetBinContent(r_bin));
     		r_bin_center_bar.push_back(hist_rBbar.GetBinCenter(r_bin));
     		r_bin_cont_bar.push_back(hist_rBbar.GetBinContent(r_bin));
   	}


	TFile fres("results_tuple_plus.root", "readonly");
	TFile fresbar("results_tuple_minus.root", "readonly");

	TTree *treeres = dynamic_cast<TTree*>(fres.Get("restree1"));
	TTree *treeresbar = dynamic_cast<TTree*>(fresbar.Get("restree1"));

	int nall = (7+7)*6 +6; 
	double res[nall], resbar[nall];
	for (int ivar=0; ivar<nall; ivar++){
		treeres->SetBranchAddress(Form("res%d", ivar), &res[ivar]);
		treeresbar->SetBranchAddress(Form("res%d", ivar), &resbar[ivar]);
	}

	treeres->GetEntry(Config::randSeed);	
	treeresbar->GetEntry(Config::randSeed);	

	sWeights sW(Config::input_file.c_str());
	auto vec_sWeights = sW.get_sWeigths(res, true);

	auto vec_sWeights_bar = sW.get_sWeigths(resbar, false);
	auto sum_sWeights = vec_sWeights;
    	// Append elements of vec2 to result
    	sum_sWeights.insert(sum_sWeights.end(), vec_sWeights_bar.begin(), vec_sWeights_bar.end());
	TRandom3 rand(Config::randSeed);
	int emax = (int)sum_sWeights.size();
        for (int e=0; e<emax; e++){ 			//Blinding
		double draw=rand.Uniform(0.0, 1.0);
		if (draw<0.5){	
			hist_B->Fill(std::get<0>(sum_sWeights[e]), std::get<1>(sum_sWeights[e]));
			logHistB->Fill(std::get<0>(sum_sWeights[e]), std::get<1>(sum_sWeights[e]));
		}else{
			hist_Bbar->Fill(std::get<0>(sum_sWeights[e]), std::get<1>(sum_sWeights[e]));
			logHistBbar->Fill(std::get<0>(sum_sWeights[e]), std::get<1>(sum_sWeights[e]));		    }	

	}
	/*TFile fsweights("draw_hist_sweighted.root", "readonly");
	//TFile fsweights("draw_hist_sweighted07.root", "readonly");
        TH1D *hist_B = dynamic_cast<TH1D *>(fsweights.Get("hist_B"));
        TH1D *hist_Bbar = dynamic_cast<TH1D *>(fsweights.Get("hist_Bbar"));
        TH1D *logHistB = dynamic_cast<TH1D *>(fsweights.Get("logHistB"));
        TH1D *logHistBbar = dynamic_cast<TH1D *>(fsweights.Get("logHistBbar"));
*/

 

        auto func_conv = [&PDFs, r_bin_center, r_bin_cont, Nrbins, Apol, Bpol, Cpol] (double x, const double *par, bool bbar) -> double{
  		double value =0.0;
   		double Tau = x;
   		for (int r_bin=1; r_bin<=Nrbins; ++r_bin){
                        double r_rec;
                        double n_r;
                        if(bbar){
                                r_rec = r_bin_center[r_bin-1];
                                n_r = r_bin_cont[r_bin-1];
                        }else{

                                r_rec = r_bin_center[r_bin-1];
                                n_r = r_bin_cont[r_bin-1];
                        }


                        double time_c =r_rec*(Apol*Tau*Tau+Bpol*Tau+Cpol); //time of mean exp value in bin     			
     			if (time_c<Params::tMin_gen || time_c> Params::tMax_gen ) continue;
     			if (bbar) 
				value += PDFs.time_acceptance(&time_c, &par[4])*PDFs.func1D_f(&time_c, par)*n_r;
			else
				value += PDFs.time_acceptance(&time_c, &par[7])*PDFs.func1D_fbar(&time_c, par)*n_r;
		}
   	return value;
	};
	auto func_norm = [func_conv] (double *x, double *par)->double{ 
		return func_conv(x[0], par, true);
	};
	auto func_norm_bar = [func_conv] (double *x, double *par)->double{ 
		return func_conv(x[0], par, false);
	};

	gRandom->SetSeed(Config::randSeed);
	double imz_gen = Config::imz_start;

	std::cout << " IIIZZZZMMMMMM " << imz_gen << std::endl << std::endl;
	//double params[10] = {0.0, imz_gen, 0.657552, 0.997031, -0.250837, 1.10646, 0.056071, -0.250837, 1.10646, 0.056071};
	TF1 f_norm("f_norm", func_norm, Params::tMin, Params::tMax, 10);
	TF1 f_norm_bar("f_norm_bar", func_norm_bar, Params::tMin, Params::tMax, 10);
	
		TH1D  *logAsym;  

			TH1D* logAsymtemp = dynamic_cast<TH1D *>(logHistB->GetAsymmetry(logHistBbar));
			logAsym = dynamic_cast<TH1D *>(logAsymtemp->Clone("logAsym"));
			TH1D *logExpAsym = dynamic_cast<TH1D *>(logAsym->Clone("logExpAsym"));
			logExpAsym->Reset();
	std::cout << hist_B->Integral() << "  " << hist_Bbar->Integral() << " Integrals signal\n\n";
	filetree.cd();
	hist_B->Write();
	hist_Bbar->Write();
        

	//hist_expect[0] = dynamic_cast<TH1D*>(hist_B->Clone("hExpect0"));
	//hist_expect[1] =  dynamic_cast<TH1D*>(hist_Bbar->Clone("hExpect1"));
	//auto logHist_expectB = dynamic_cast<TH1D *>(logHistB->Clone("expectB"));
	//auto logHist_expectBbar = dynamic_cast<TH1D *>(logHistBbar->Clone("expectB"));

	filetree.cd();
	double output[10], chi2, error[10];
	
	TTree *outtree[2];
        outtree[0] = new TTree("outtree1", "outtree1");
        outtree[1] = new TTree("outtree2", "outtree2");

	for (int ii=0; ii<2; ii++){
	outtree[ii]->Branch("imz_gen", &imz_gen, "imz_gen/D");
	outtree[ii]->Branch("chi2", &chi2, "chi2/D");
	outtree[ii]->Branch("rez", &output[0], "rez/D");
	outtree[ii]->Branch("drez", &error[0], "drez/D");
	outtree[ii]->Branch("imz", &output[1], "imz/D");
	outtree[ii]->Branch("dimz", &error[1], "dimz/D");
	outtree[ii]->Branch("Gamma", &output[2], "Gamma/D");
	outtree[ii]->Branch("dGamma", &error[2], "dGamma/D");
	outtree[ii]->Branch("eta", &output[3], "eta/D");
	outtree[ii]->Branch("deta", &error[3], "deta/D");
	outtree[ii]->Branch("t_shift", &output[4], "t_shift/D");
	outtree[ii]->Branch("dt_shift", &error[4], "dt_shift/D");
	outtree[ii]->Branch("alpha", &output[5], "alpha/D");
	outtree[ii]->Branch("dalpha", &error[5], "dalpha/D");
	outtree[ii]->Branch("beta", &output[6], "beta/D");
	outtree[ii]->Branch("dbeta", &error[6], "dbeta/D");
	outtree[ii]->Branch("t_shift_anit", &output[7], "t_shift_anti/D");
	outtree[ii]->Branch("dt_shift_anit", &error[7], "dt_shift_anti/D");
	outtree[ii]->Branch("alpha_anti", &output[8], "alpha_anti/D");
	outtree[ii]->Branch("dalpha_anti", &error[8], "dalpha_anti/D");
	outtree[ii]->Branch("beta_anti", &output[9], "beta_anti/D");
	outtree[ii]->Branch("dbeta_anti", &error[9], "dbeta_anti/D");

	}



	TRandom3 randvar(Config::randSeed);
	double minchi2 = 1.0e14;
        double minoutput[Config::nvar_time];
        for (int i=0; i<Config::nvar_time; i++) minoutput[i] = -1000.0;
        for (int ii=0; ii<2; ii++){
		bool make_asym = ii==1;
        for (int itry=0; itry<Config::ntries[ii]; itry++){


	filetree.cd();

        std::string minName = "Minuit2";
        std::string algoName = "";
	ROOT::Math::Minimizer *min =
                        ROOT::Math::Factory::CreateMinimizer(minName, algoName);

        min->SetMaxFunctionCalls(Config::functionCalls); // for Minuit/Minuit2
        min->SetTolerance(Config::tolerance[ii]);
        min->SetPrintLevel(Config::printLevel);
	double CL_normal = ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1); //1 sigma ~68%
	min->SetErrorDef(TMath::ChisquareQuantile(CL_normal, 4));

	
	int ivar=0;
	for (const auto& pair : Config::varname){ 
		double drawvar=1.0;
		if (ivar>=3 && ivar <=6)
			drawvar = randvar.Uniform(Config::draw_min, Config::draw_max);	
		min->SetVariable(ivar, pair.first.c_str(), drawvar*(pair.second), 0.0001);
		ivar++;
	}

	if (ii==1){
		for(int i=0; i<Config::nvar_time; i++)
			min->SetVariableValue(i, minoutput[i]);
	}
	//for (int i=0; i<10; i++) min->FixVariable(i);
	std::vector<std::pair<int, int>> replaceIndexVect = {};

        const auto& replace_var = Config::replace_var;
                for (const auto& rep_var: replace_var){
                        int index_replaced = min->VariableIndex(rep_var.first);
                        int index_substitute = min->VariableIndex(rep_var.second);
                        if (index_replaced == -1 ){
                                std::cerr<< "Error in substituting parameters param " << rep_var.first << " not found.\n";
                                return 1;
                        }
                        if (index_substitute == -1 ){
                                std::cerr<< "Error in substituting parameters param " << rep_var.second << " not found.\n";
                                return 1;
                        }
			Config::fixVect[ii].push_back(rep_var.first);
        	        replaceIndexVect.push_back(std::make_pair(index_replaced, index_substitute));
                }

	for (const auto& fix: Config::fixVect[ii]){
                min->FixVariable(min->VariableIndex(fix));
        }



	auto func_chi2 = wrap_chi2(logHist_expectB, logHist_expectBbar, logAsym, logExpAsym, func_conv, vec_time, vec_time_bar, hist_B, hist_Bbar, logHistB, logHistBbar, make_asym, replaceIndexVect);	

	ROOT::Math::Functor f(func_chi2, Config::nvar_time);
	min->SetFunction(f);
        min->Minimize();

	chi2 = min->MinValue();

	filetree.cd();
	const double *out = min->X();
	for (int i=0; i<Config::nvar_time; i++)
		output[i] = out[i];
	const double *err = min->Errors();
	for (int i=0; i<Config::nvar_time; i++)
		error[i] = err[i];
	        
	outtree[ii]->Fill();
        if (ii==0 && minchi2>chi2){
        	for (int i=0; i<Config::nvar_time; i++) minoutput[i] = output[i];
      	        minchi2 = chi2;
        }
 
		

	       		       
	if (min) delete min;
}
}
	outtree[0]->Write();
	outtree[1]->Write();


	filetree.cd();
	for (int int_bbar=1; int_bbar<2; int_bbar++){
			for (int bin=1; bin<=Config::nlogBins; bin++){
				double t = logHist_expectB->GetBinCenter(bin);
				double expect = func_conv(t, output, int_bbar==0);
				if (int_bbar==0){
					logHist_expectB->SetBinContent(bin, expect*logHist_expectB->GetBinWidth(bin));
					logHist_expectB->SetBinError(bin, 0.0);
				}else{
					logHist_expectBbar->SetBinContent(bin, expect*logHist_expectBbar->GetBinWidth(bin));
					logHist_expectBbar->SetBinError(bin, 0.0);
				}
			}
		}
		logHist_expectB->Scale(logHistB->Integral()/logHist_expectB->Integral());
		logHist_expectBbar->Scale(logHistBbar->Integral()/logHist_expectBbar->Integral());
	//}
	logHistB->Write();
	logHistBbar->Write();
	logHist_expectB->Write();
	logHist_expectBbar->Write();


	TH1D *logExpAsymtemp =  dynamic_cast<TH1D *>(logHist_expectB->GetAsymmetry(logHist_expectBbar));
	logExpAsym = dynamic_cast<TH1D *>(logExpAsymtemp->Clone("logAsym_res"));
		logExpAsym->Write();
		logAsym->Write();
	
				
	filetree.Close();

	return 0;
}
std::function<double(const double*)> wrap_chi2(TH1D *logHist_expectB, TH1D *logHist_expectBbar, TH1D *logAsym, TH1D *logExpAsym, std::function<double(double, const double*, bool)> func_conv, std::vector<double> vec_time, std::vector<double> vec_time_bar,  TH1D *hist_B, TH1D *hist_Bbar, TH1D *logHistB, TH1D *logHistBbar, bool make_asym, std::vector<std::pair<int, int>> replaceIndexVect){
	auto func_chi2 = [logHist_expectB, logHist_expectBbar, logAsym, func_conv, vec_time, vec_time_bar, hist_B, hist_Bbar, logHistB, logHistBbar, make_asym, replaceIndexVect] (const double *par)->double{
		// Extract the parameters and add some constraints
                double param[Config::nvar_time];
                for (int ivar = 0; ivar < Config::nvar_time; ivar++)
                {
                        param[ivar] = par[ivar];
                }
                for (const auto& irep_var: replaceIndexVect){
                        param[irep_var.first] = par[irep_var.second];
                }





		double chi2 = 0.0;

		for (int int_bbar=0; int_bbar<2; int_bbar++){
				double norm = 0.0;
				std::vector<double> vect_expect;
				vect_expect.clear();
				for (int bin=1; bin<=Config::nlogBins; bin++){
                                        double t = logHistB->GetBinCenter(bin);

                                        double expect = func_conv(t, param, int_bbar==0)*logHistB->GetBinWidth(bin);
					vect_expect.push_back(expect);
					norm +=expect;
                                }

				double integral = hist_B->Integral();
				double integralBbar = hist_Bbar->Integral();
				for (int bin=1; bin<=Config::nlogBins; bin++){
					if(int_bbar==1) 
						integral = integralBbar;
					double expect = vect_expect[bin-1]/norm*integral;
					if(!make_asym){
						double cont, err;
						if (int_bbar==0){
							cont = logHistB->GetBinContent(bin); 
							err = logHistB->GetBinError(bin);
						}else{	
							cont = logHistBbar->GetBinContent(bin); 
							err = logHistBbar->GetBinError(bin);
						}
						if (err!=0)
							chi2+= (cont-expect)*(cont-expect)/err/err;
					}
					if (int_bbar==0){
						logHist_expectB->SetBinContent(bin, expect);
						logHist_expectB->SetBinError(bin, 0.0);
					}else{
						logHist_expectBbar->SetBinContent(bin, expect);
						logHist_expectBbar->SetBinError(bin, 0.0);
					}
					
				}
			
				if (make_asym){
					if(int_bbar==0) continue;
					//logHist_expectB->Scale(logHistB->Integral()/logHist_expectB->Integral());
					//logHist_expectBbar->Scale(logHistBbar->Integral()/logHist_expectBbar->Integral());
					TH1D *logExpAsym = dynamic_cast<TH1D*>(logHist_expectB->GetAsymmetry(logHist_expectBbar)); 	
					for (int bin=1; bin<=Config::nlogBins; bin++){
						double cont = logAsym->GetBinContent(bin);
						double err = logAsym->GetBinError(bin);
						double expect = logExpAsym->GetBinContent(bin);
						if (err!=0)
							chi2+= (cont-expect)*(cont-expect)/err/err;
					}
				}
		}

		
		return chi2;
	};
	return func_chi2;
}
