#include <fstream>
#include <iomanip>
#include <string>

#include "config.h"
#include "json.hpp"
#include "omp.h"
#include "sweights.h"

// ROOT includes
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/ProbFuncMathCore.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TString.h"

// CPT_beauty_2Dmass_fit
#include <functional>
#include <memory>
#include <tuple>

#include "ChebyshevPDF.h"
#include "D_M_fit_shape.h"
#include "FastSum.h"
#include "M_B_2missPT_fit.h"
#include "TH2D.h"
#include "TSystem.h"

using json = nlohmann::json;
using namespace cpt_b0_analysis;

#include <mutex>
std::mutex my_mutex;

// bool avx = false;
bool avx = true;

std::function<double(const double*)> wrap_chi2(
	const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get,
	const std::vector<std::shared_ptr<PDFInterface>>& B_PDFs_get,
	std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD,
	const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB,
	const std::vector<std::vector<double>>& dMC_MB,
	const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH2D hist2D);

std::function<double(const double*)> wrap_chi2_simultanous(
	const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get,
	const std::vector<std::shared_ptr<PDFInterface>>& B_PDFs_get,
	std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD,
	const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB,
	const std::vector<std::vector<double>>& dMC_MB,
	const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH1D hist_MD, TH1D hist_MB);

/**
 * @brief A 2D fitter for B and D mass distributions
 *
 * @param argc
 * @param argv - takes a single argument - path to a config file
 * @return int
 */
int main(int argc, char* argv[]) {
	// The first, fundamental operation to be performed in order to make ROOT
	// thread-aware.
	ROOT::EnableThreadSafety();

	// Increase printout precision
	// std::cout<<std::setprecision(40);

	// Load config
	try {
		if (argc != 2) {
			std::cerr << "Usage: fit2D_mass config.json" << std::endl;
			return 1;
		}
		if (Config::load(argv[1])) {
			std::cerr << " Bad config file! " << std::endl;
			return 1;
		}

	} catch (const std::exception& ex) {
		std::cerr << "Error: " << ex.what() << std::endl;
		return 1;
	}
	// Limit number of events to read
	int nentries = Config::nentries;

	// Load data set
	TChain ch((Config::chainName).c_str());
	ch.Add(Config::input_files[0].c_str());
	double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT, Tau;
	double B_MMcorr;
	bool charge;
	int frac_index = 100;  // TODO_DOCS - what is this ??
	const int nbins = 100;
	// TODO_KK 2x1D version
	TH2D hist2D("hist2D", "", nbins, Config::minDM, Config::maxDM, nbins, Config::minBMcorr,
				Config::maxBMcorr);
	TH1D hist_MD("hist_MD", "", nbins, Config::minDM, Config::maxDM);
	TH1D hist_MB("hist_MB", "", nbins, Config::minBMcorr, Config::maxBMcorr);

	std::vector<std::pair<double, double>> vect_2D = {};
	if (!Config::isMC) {
		ch.SetBranchAddress("B_M", &B_M);
		ch.SetBranchAddress("missPT", &missPT);
		ch.SetBranchAddress("D_M", &D_M);
		ch.SetBranchAddress("mu_PT", &mu_PT);
		ch.SetBranchAddress("mu_P", &mu_P);
		ch.SetBranchAddress("mu_eta", &mu_eta);
		ch.SetBranchAddress("K_PT", &K_PT);
		ch.SetBranchAddress("Tau", &Tau);
		ch.SetBranchAddress("truecharge", &charge);
	} else {
		ch.SetBranchAddress("BMcorr", &B_MMcorr);
		ch.SetBranchAddress("DM", &D_M);
		ch.SetBranchAddress("frac_index", &frac_index);
	}

	if (nentries < 0) nentries = ch.GetEntries();
	if (nentries > ch.GetEntries()) {
		std::cerr << "The value of 'nentries' exceeds the number of events in the file."
				  << std::endl;
		return 1;
	}
	int frac_indeces[Config::ncontr];
	for (int i = 0; i < Config::ncontr; i++) frac_indeces[i] = 0;
	std::cout << nentries << " number of entries read" << std::endl;
	for (int i = 0; i < nentries; ++i) {
		ch.GetEntry(i);
		// Apply cuts
		if (!Config::isMC) {
			if (mu_PT < Config::muPTmin || mu_P < Config::muPmin || mu_eta < Config::eta_min ||
				mu_eta > Config::eta_max)
				continue;
			if (int(charge) != Config::sign) continue;
			B_MMcorr = B_M + 2.0 * missPT;
		}
		if (B_MMcorr < Config::minBMcorr || B_MMcorr > Config::maxBMcorr) continue;

		if (D_M < Config::minDM || D_M > Config::maxDM) continue;

		if (Tau < Config::tMin || Tau > Config::tMax) continue;

		vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
		hist2D.Fill(D_M, B_MMcorr);
		hist_MD.Fill(D_M);
		hist_MB.Fill(B_MMcorr);

		// TODO_DOCS count number of entries for given contribution ????
		if (Config::isMC) frac_indeces[frac_index]++;
	}
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cout << hist2D.Integral() << "  " << vect_2D.size() << std::endl;

	// Initial fit parameter values taken from 1D fits to MC and Side Bands

	// Define Minuit fit variables for M_D
	const int n_all = Config::nvar_all_md + Config::nvar_all_mb + Config::ncontr;
	// Get index of the sidebands contribution
	const int sb_idx =
		std::distance(Config::contrName.begin(),
					  std::find(Config::contrName.begin(), Config::contrName.end(), "sidebands"));
	// Indicates that the previous try was a good fit
	bool previous_fit = false;
	double starting_point[n_all];
	std::string minName = "Minuit2";
	std::string algoName = "";
	int itry = 1;
	// Indicates that at least one fit converged
	bool goodfit = false;
	// Currently selected fit ID from "Fits" config attibute
	int last_fit = -1;
	// Repeat fits until a good one is achieved
	while ((itry <= Config::ntries || !goodfit) && itry <= 100) {
		bool start_scratch = true;
		int loop_count = 0;
		for (auto& int_choose_fit : Config::int_choose_fits) {
			last_fit = int_choose_fit;
			ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

			// Set tolerance , etc...
			min->SetMaxFunctionCalls(Config::functionCalls);  // for Minuit/Minuit2
			min->SetTolerance(Config::tolerance[loop_count]);
			loop_count++;
			min->SetPrintLevel(Config::printLevel);
			// min->SetStrategy(2);

			// Define M_D PDF parameters
			for (int i = 0; i < Config::ncontr; i++) {
				for (int ivar = 0; ivar < Config::nvar_md[i]; ivar++) {
					// Set initial values from 1D fit best value
					// Set step size based on the 1D fit uncertainty
					min->SetVariable(Config::nvar_offset_md[i] + ivar,
									 (TString::Format("md_%s_%s", Config::contrName[i].c_str(),
													  Config::varname_md[i][ivar].c_str()))
										 .Data(),
									 Config::MC_MD[i][ivar], Config::dMC_MD[i][ivar] + 1.0e-11);
					// For BM and frac fits fix all D_M PDF params
					if (((int_choose_fit == dictionaryChooseFit.at("BM") ||
						  int_choose_fit == dictionaryChooseFit.at("frac")) &&
						 i != -1)) {
						min->FixVariable(Config::nvar_offset_md[i] + ivar);
					}
					// Store initial values for "from scratch" fit
					if (start_scratch)
						starting_point[Config::nvar_offset_md[i] + ivar] = Config::MC_MD[i][ivar];
				}
			}

			// Define Minuit fit variables for M_B

			for (int i = 0; i < Config::ncontr; i++) {
				for (int ivar = 0; ivar < Config::nvar_mb[i]; ivar++) {
					// Set initial values from 1D fit best value
					// Set step size based on the 1D fit uncertainty
					min->SetVariable(Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar,
									 (TString::Format("mb_%s_%s", Config::contrName[i].c_str(),
													  Config::varname_mb[i][ivar].c_str()))
										 .Data(),
									 Config::MC_MB[i][ivar], Config::dMC_MB[i][ivar] + 1.0e-11);
					// For DM+BMfixed and frac fits fix all B_M PDF parameters
					if (int_choose_fit == dictionaryChooseFit.at("DM+BMfixed") ||
						int_choose_fit == dictionaryChooseFit.at("frac")) {
						min->FixVariable(Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar);
					}
					// Store initial values for "from scratch" fit
					if (start_scratch)
						starting_point[Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar] =
							Config::MC_MB[i][ivar];
				}
			}

			// Set Limits on variables from the config file
			for (auto it = Config::varLimitsMap.begin(); it != Config::varLimitsMap.end(); ++it) {
				int index_var = min->VariableIndex(it->first);
				if (index_var == -1) {
					std::cerr << "Error in limiting parameters: param " << it->first
							  << " not found.\n";
					return 1;
				}
				auto pairlims = it->second;
				min->SetVariableLimits(index_var, pairlims.first, pairlims.second);
			}

			// Load sideband fraction from fit results
			double frac_sidebands = abs(Config::MC_MD[sb_idx][Config::nvar_md[sb_idx]]);
			Config::fracInit[sb_idx] = frac_sidebands;
			// Define contribution fractions
			for (int i = 0; i < Config::ncontr; i++) {
				// Load initial fractions from config
				min->SetVariable(Config::nvar_all_md + Config::nvar_all_mb + i,
								 (TString::Format("par_frac%d", i)).Data(), Config::fracInit[i],
								 0.001);
				if (start_scratch)
					starting_point[Config::nvar_all_md + Config::nvar_all_mb + i] =
						Config::fracInit[i];
			}
			// Fix the last fraction. It will be calculated on the fly from normalisation to 1.
			min->FixVariable(Config::nvar_all_md + Config::nvar_all_mb + Config::ncontr - 1);

			// Define the error setimation parameter in minuit for 1 sigma and ncontr -1 free
			// parameters
			// TODO_KK - make it configurable via config
			double CL_normal =
				ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1);	 // 1 sigma ~68%
			min->SetErrorDef(TMath::ChisquareQuantile(
				CL_normal,
				Config::ncontr - 1));  // ncontr-1 free fraction parameters, other parameters have
									   // gaussian contraints base on MC fits.
			
						

			// Load D_M PDFs
			const auto& D_PDFs = Config::getVectorPDFs("Dmass");
			if (int(D_PDFs.size()) != Config::ncontr) {
				std::cout << " NO D_PDFs \n";
				return 1;
			}
			for (int i = 0; i < Config::ncontr; ++i) {
				auto pdf = D_PDFs[i].get();

				if (!pdf) {
					std::cout << "Nullptr passed as pdf\n";
					return 1;
				}
			}

			// Load B_M PDFs
			const auto& B_PDFs = Config::getVectorPDFs("Bmass");
			if (int(B_PDFs.size()) != Config::ncontr) {
				std::cout << " NO B_PDFs \n";
				return 1;
			}
			for (int i = 0; i < Config::ncontr; ++i) {
				auto pdf = B_PDFs[i].get();
				if (!pdf) {
					std::cout << "Nullptr passed as pdf\n";
					return 1;
				}
			}

			// Fix variables form config
			for (const auto& fix : Config::fixVect) {
				min->FixVariable(min->VariableIndex(fix));
				// std::cout << fix << "  " << min->VariableIndex(fix) << std::endl;
			}

			// Set random seed
			TRandom rand;
			if (Config::randSeed > -1) rand.SetSeed(itry);

			// Randomly smear staring point for stability checks
			// Only for fit "all" with defined randSeed and not the first fit
			double frac_st[Config::ncontr];
			for (int ivar = 0; ivar < Config::nvar_all_md + Config::nvar_all_mb + Config::ncontr;
				 ivar++) {
				double random = 0.0;
				if (Config::randSeed > -1 && !min->IsFixedVariable(ivar) &&
					int_choose_fit == dictionaryChooseFit.at("all") && itry != 1)
					random = rand.Uniform(-1.0, 1.0);

				min->SetVariableValue(ivar, starting_point[ivar] * (1.0 + 0.01 * random));
			}

			// TODO_DOCS why this is hardcoded and only read from files for "frac" and "all" ???
			// TODO why not use teh Config::MD_MC ???
			// For other fits this is not set at all ???

			// TODO this is not needed if we have supprort to per component PDFs
			/* KK
			double frac_sidebands = -0.270922;
			if (int_choose_fit == dictionaryChooseFit.at("frac") ||
				int_choose_fit == dictionaryChooseFit.at("all")) {
				std::ifstream sideband_input(Form("results1D_DM_%d_1.txt", Config::sign));
				double a1 = -0.00102346;
				double a2 = 1.32508e-07;
				double dD = 0.0;
				sideband_input >> dD >> dD; // Fit status
				sideband_input >> a1 >> dD;
				sideband_input >> a2 >> dD;
				sideband_input >> frac_sidebands >> dD;
				std::cout << "========================================\n read data " << a1 << "  "
						  << a2 << "  " << frac_sidebands << std::endl;
				double a_pair[] = {a1, a2};
				for (int ivar = 0; ivar < Config::n_sideband; ivar++)
					min->SetVariableValue(Config::nvar_offset_md[sb_idx] + ivar, a_pair[ivar]);
				min->SetVariableValue(Config::nvar_all_md + Config::nvar_all_mb + sb_idx,
									  frac_sidebands);
			}
			*/

			double sumc = frac_sidebands;
			// TODO rewrite add randomFix config param
			for (int icontr = 0; icontr < Config::ncontr; icontr++) {
				frac_st[icontr] =
					starting_point[Config::nvar_all_md + Config::nvar_all_mb + icontr];

				if (icontr == sb_idx) continue;

				if (Config::randSeed > -1 &&
					!min->IsFixedVariable(Config::nvar_all_md + Config::nvar_all_mb + icontr) &&
					start_scratch && !Config::start_from_previous && itry != 1) {
					if (icontr == Config::ncontr - 1) {
						frac_st[icontr] = 1 - sumc;
						continue;
					}
					double random = rand.Uniform(0.0, 1.0 - sumc);
					frac_st[icontr] = random;

					sumc += abs(frac_st[icontr]);
				}
			}

			double result0[Config::nvar_all_md + Config::nvar_all_mb + 2 * Config::ncontr];
			if (Config::start_from_previous) {
				double x, dx;
				std::ifstream res0(Config::previous_result_file.c_str());
				// Read fit status and chi2/NLL value
				res0 >> x >> dx;
				int count = 0;
				while (res0 >> x >> dx) {
					result0[count] = x;
					count++;
				}
				for (int ic = 0; ic < Config::ncontr; ic++) {
					double random = 0.0;
					if (Config::randSeed > -1 &&
						!min->IsFixedVariable(Config::nvar_all_md + Config::nvar_all_mb + ic) &&
						(int_choose_fit == dictionaryChooseFit.at("all") ||
						 int_choose_fit == dictionaryChooseFit.at("frac")) &&
						itry != 1)
						random = rand.Uniform(-1.0, 1.0);
					frac_st[ic] = result0[Config::nvar_all_md + Config::nvar_all_mb + ic] *
								  (1.0 + 0.01 * random);
				}
			}
			for (int icontr = 0; icontr < Config::ncontr; icontr++)
				std::cout << frac_st[icontr] << "  <===frac" << icontr << std::endl;

			// Fix small contributions
			double sum_contr = 0.0;
			for (int ic = 0; ic < Config::ncontr; ic++) {
				sum_contr += abs(frac_st[ic]);
				if (abs(frac_st[ic]) < 1.0e-4) {
					std::cout << ic << "  fixed\n";
					// frac_st[ic] = 0.0;
					for (int ivar = 0; ivar < Config::nvar_md[ic]; ivar++)
						min->FixVariable(Config::nvar_offset_md[ic] + ivar);
					for (int ivar = 0; ivar < Config::nvar_mb[ic]; ivar++)
						min->FixVariable(Config::nvar_all_md + Config::nvar_offset_mb[ic] + ivar);
					// min->FixVariable((Config::nvar_md+Config::nvar_mb) * Config::ncontr+ic);
				}
			}

			// Make sure that fractions are normalized ...
			for (int ic = 0; ic < Config::ncontr; ic++) {
				min->SetVariableValue(Config::nvar_all_md + Config::nvar_all_mb + ic,
									  frac_st[ic] / sum_contr);
				// Fix fractions for shape fits
				if (int_choose_fit == dictionaryChooseFit.at("shapes") ||
					int_choose_fit == dictionaryChooseFit.at("DM+BMfixed") ||
					int_choose_fit == dictionaryChooseFit.at("BM"))
					min->FixVariable(Config::nvar_all_md + Config::nvar_all_mb + ic);
			}

			// Generate variable substitution rules
			std::vector<std::pair<int, int>> replaceIndexVect = {};
			for (const auto& rep_var : Config::replace_var) {
				int index_replaced = min->VariableIndex(rep_var.first);
				int index_substitute = min->VariableIndex(rep_var.second);
				if (index_replaced == -1) {
					std::cerr << "Error in substituting parameters param " << rep_var.first
							  << " not found.\n";
					return 1;
				}
				if (index_substitute == -1) {
					std::cerr << "Error in substituting parameters param " << rep_var.second
							  << " not found.\n";
					return 1;
				}
				replaceIndexVect.push_back(std::make_pair(index_replaced, index_substitute));
			}

			// Start the minimization
			// Define a fit function for Minuit
			bool simulatnous = true;
			std::function<double(const double*)> fchi2;
			if (simulatnous) {
				fchi2 =
					wrap_chi2_simultanous(D_PDFs, B_PDFs, vect_2D, Config::MC_MD, Config::dMC_MD, Config::MC_MB,
							Config::dMC_MB, replaceIndexVect, int_choose_fit, hist_MD, hist_MB);
			} else {
				fchi2 =
					wrap_chi2(D_PDFs, B_PDFs, vect_2D, Config::MC_MD, Config::dMC_MD, Config::MC_MB,
							  Config::dMC_MB, replaceIndexVect, int_choose_fit, hist2D);
			}
			ROOT::Math::Functor f(fchi2, n_all);
			// ROOT::Math::Functor f(fchi2, (nvar_md + nvar_mb) * ncontr + ncontr);
			min->SetFunction(f);
			min->Minimize();

			// TODO get output path from config
			// Print the fit results
			TString path, path1;
			if (Config::binned) {
				path1 = "results_binned100kMU_2";
				path = TString::Format("results_binned100kMU_2/fit2D_%d", itry);
			} else {
				path1 = "results_unbinned100kMU_2";
				path = TString::Format("results_unbinned100kMU_2/fit2D_%d", itry);
			}
			gSystem->Exec(TString::Format("mkdir -p %s", path1.Data()).Data());
			gSystem->Exec(TString::Format("mkdir -p %s", path.Data()).Data());
			std::ofstream results(
				TString::Format("%s/results_%d_%d.txt", path.Data(), Config::sign, int_choose_fit));
			std::cout << std::setprecision(25);
			results << std::setprecision(25) << min->Status() << "  " << min->MinValue()
					<< std::endl;
			std::cout << std::setprecision(10);
			if (min->Status() != 0 && min->Status() != 1) {	 // || min->MinValue()>24.8e6){
				std::cout << "Bad status of fit " << min->Status() << " chi2 is " << min->MinValue()
						  << std::endl;
				previous_fit = false;
			} else {
				goodfit = true;
				previous_fit = true;
				start_scratch = false;
			}
			for (int i = 0; i < n_all - 1; i++) {
				results << min->X()[i] << "  " << min->Errors()[i] << std::endl;
				if (previous_fit) starting_point[i] = min->X()[i];
			}
			const double* pa = &min->X()[Config::nvar_all_md + Config::nvar_all_mb];
			double frac_res[Config::ncontr];
			double sumfrac = 0.0;
			for (int i = 0; i < Config::ncontr - 1; i++) {
				frac_res[i] = abs(pa[i]);
				sumfrac += frac_res[i];
			}
			double frac_last = 1.0 - sumfrac;
			frac_res[Config::ncontr - 1] = frac_last;
			results << frac_last << "  " << 0.0 << std::endl;
			if (previous_fit)
				starting_point[Config::nvar_all_md + Config::nvar_all_mb + Config::ncontr - 1] =
					frac_last;

			std::string fileWeightsName = "Tree_sWeights";
			std::string TreeName = "Tree_sWeights";
			// TODO_KK fix sWeights
			//sWeights sW(Config::input_files[0].c_str(), fileWeightsName.c_str(), TreeName.c_str());
			//sW.get_sWeigths(min->X(), Config::sign == 1);

			// TODO_DOCS What is this ?? Why are we storing this value (effectively read from flags from ROOT) and not the fit results? 
			for (int i = 0; i < Config::ncontr; i++)
				results << double(frac_indeces[i]) / double(vect_2D.size()) << "  " << 0.0
						<< std::endl;

			results.close();

			if (min) {
				delete min;
			}
			for (int i = 0; i < Config::ncontr; i++)
				std::cout << frac_res[i] << "  frac" << i << "  "
						  << double(frac_indeces[i]) / double(vect_2D.size()) << std::endl;
		}
		itry++;
	}
	double min_chi2 = 1.0e45;
	double chi2;
	int status;
	int best = -1;
	TString result_dir;
	if (Config::binned)
		result_dir = "results_binned100kMU_2";
	else
		result_dir = "results_unbinned100kMU_2";
	for (int i = 1; i <= itry; i++) {
		double res[Config::nvar_all_md + Config::nvar_all_mb + 2 * Config::ncontr];
		double dres[Config::nvar_all_md + Config::nvar_all_mb + 2 * Config::ncontr];
		std::ifstream resin(
			Form("%s/fit2D_%d/results_%d_%d.txt", result_dir.Data(), i, Config::sign, last_fit));
		resin >> status >> chi2;

		int j = 0;
		double x, dx;
		while (resin >> x >> dx) {
			res[j] = abs(x);
			dres[j] = dx;
			j++;
		}

		resin.close();

		int n0 = Config::nvar_all_md + Config::nvar_all_mb;
		if ((status == 0 || status == 1) &&
			chi2 < min_chi2) {	// && res[n0+1]<0.2 && res[n0+2] <0.2 && res[n0+3]<0.2&&
								// res[n0+5]<0.2){// && i!=10){
			min_chi2 = chi2;
			best = i;
		}
	}
	TString path, best_path;
	if (Config::binned) {
		path = TString::Format("results_binned100kMU_2/fit2D_%d/", best);
		best_path = "best_results_binned/fit2D_best";
	} else {
		path = TString::Format("results_unbinned100kMU_2/fit2D_%d/", best);
		best_path = "best_results_unbinned/fit2D_best/";
	}
	gSystem->Exec(TString::Format("mkdir -p  %s", best_path.Data()).Data());
	gSystem->Exec(TString::Format("cp -r  %s/* %s", path.Data(), best_path.Data()).Data());

	return 0;
}

std::function<double(const double*)> wrap_chi2(
	const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get,
	const std::vector<std::shared_ptr<PDFInterface>>& B_PDFs_get,
	std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD,
	const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB,
	const std::vector<std::vector<double>>& dMC_MB,
	const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH2D hist2D) {
	long int emax = vect_2D.size();
	std::cout << emax << " yield of the sample \n";
	auto fchi2 = [D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB, emax,
				  replaceIndexVect, int_choose_fit, hist2D](const double* par) -> double {
		// Caclulate chi2
		// Get the pointer in parameters array that correspond to fraction defining parameters????
		const double* pa = &par[Config::nvar_all_md + Config::nvar_all_mb];
		// Calculate fractions
		// The parameters pa[] define the fractions, we have 6 fractions but 5 independent
		// parameters. The parametrisation is arbitrary
		double frac[Config::ncontr];

		double chi2 = 0.0;
		double sum_frac = 0.0;
		for (int i = 0; i < Config::ncontr - 1; i++) {
			frac[i] = abs(pa[i]);
			sum_frac += frac[i];
		}
		frac[Config::ncontr - 1] = abs(1.0 - sum_frac);
		sum_frac += frac[Config::ncontr - 1];
		double tmp = 1.0e10 * (sum_frac - 1.0) * (sum_frac - 1.0);
		if (!Config::binned) tmp *= emax;
		chi2 += tmp;

		/*frac[0] = 1.0 - abs(pa[0]);
		frac[1] = abs(pa[0]) * (1.0 - abs(pa[1]));
		frac[2] = abs(pa[0]) * abs(pa[1]) * (1.0 - abs(pa[2]));
		frac[3] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * (1.0 - abs(pa[3]));
		frac[4] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * abs(pa[3]) * (1.0 - abs(pa[4]));
		frac[5] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * abs(pa[3]) * abs(pa[4]);
	*/
		// Extract the parameters and add some constraints
		double param[Config::nvar_all_md + Config::nvar_all_mb];
		for (int ivar = 0; ivar < (Config::nvar_all_md + Config::nvar_all_mb); ivar++) {
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var : replaceIndexVect) {
			param[irep_var.first] = par[irep_var.second];
		}

		// Calculate normalisation integrals
		for (int i = 0; i < Config::ncontr; i++) {
			D_PDFs_get[i]->CalcIntegral(&param[Config::nvar_offset_md[i]], Config::minDM,
										Config::maxDM);
			B_PDFs_get[i]->CalcIntegral(&param[Config::nvar_all_md + Config::nvar_offset_mb[i]],
										Config::minBMcorr, Config::maxBMcorr);
		}

		int nbins_md = hist2D.GetNbinsX();
		int nbins_mb = hist2D.GetNbinsY();
		if (int_choose_fit !=
			dictionaryChooseFit.at("frac")) {  // don't calculate if only fractions fit
			for (int i = 0; i < Config::ncontr; i++) {
				// TODO set dxx for D comb background to 0 and remove this
				for (int ivar = 0; ivar < Config::nvar_md[i]; ivar++) {
					// Skip the constrain for chebyshev background: {"Chebyshev", 1}
					// if (Config::intshapesDM[i] == 1 || ivar == 1) continue;
					if (dMC_MD[i][ivar] != 0) {
						double tmp = (MC_MD[i][ivar] - param[Config::nvar_offset_md[i] + ivar]) *
									 (MC_MD[i][ivar] - param[Config::nvar_offset_md[i] + ivar]) /
									 (2.0 * (dMC_MD[i][ivar] *
											 dMC_MD[i][ivar]));	 // use the results of MC fits
						double scale = 1.0;
						chi2 += scale * tmp;
					}
				}
			}
		}
		if (int_choose_fit !=
			dictionaryChooseFit.at("frac")) {  // don't calculate if only fractions fit
			for (int i = 0; i < Config::ncontr; i++) {
				for (int ivar = 0; ivar < Config::nvar_mb[i]; ivar++) {
					if (dMC_MB[i][ivar] != 0) {
						double tmp =
							(MC_MB[i][ivar] -
							 param[Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar]) *
							(MC_MB[i][ivar] -
							 param[Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar]) /
							(2.0 *
							 (dMC_MB[i][ivar] * dMC_MB[i][ivar]));	// use the results of MC fits

						double scale = 1.0;
						// if(int_choose_fit == dictionaryChooseFit.at("BM"))
						//	scale = 10000.0;
						chi2 += scale * tmp;
					}
				}
			}
		}

		double* vect_chi2 =
			new double[vect_2D.size()];	 // Per event results - required to efficiently calculate a
										 // Kahan compensated sum
		for (int jj = 0; jj < int(vect_2D.size()); jj++) vect_chi2[jj] = 0.0;
		if (Config::binned) {
			double bin_width_mb = (Config::maxBMcorr - Config::minBMcorr) / double(nbins_mb);
			double bin_width_md = (Config::maxDM - Config::minDM) / double(nbins_md);
			double histev = hist2D.Integral();
#pragma omp parallel for
			for (int bin_md = 1; bin_md <= nbins_md; bin_md++) {
				for (int bin_mb = 1; bin_mb <= nbins_mb; bin_mb++) {
					double sum_contr = 0.0;
					for (int i = 0; i < Config::ncontr; i++) {
						double mdass = hist2D.GetXaxis()->GetBinCenter(bin_md);
						double mcorr = hist2D.GetYaxis()->GetBinCenter(bin_mb);
						double md_val =
							D_PDFs_get[i]->EvalPDF(&mdass, &param[Config::nvar_offset_md[i]]);
						double mb_val = B_PDFs_get[i]->EvalPDF(
							&mcorr, &param[Config::nvar_all_md + Config::nvar_offset_mb[i]]);
						sum_contr +=
							histev * bin_width_md * bin_width_mb * frac[i] * mb_val * md_val;
						// std::cout << sum_contr << "  sumcontr  " << frac[i] << "  " << md_val <<
						// "  " << mb_val << std::endl;
						if (frac[i] < 0.0 || frac[i] > 1.0) {
							chi2 += 1e15;
							continue;
						}
					}
					double cont = (double)hist2D.GetBinContent(bin_md, bin_mb);
					double err = (double)hist2D.GetBinError(bin_md, bin_mb);
					if (err != 0.0) {
						vect_chi2[bin_md - 1] +=
							0.5 * (sum_contr - cont) * (sum_contr - cont) / err / err;
					}
				}
			}
		} else {

			// Main loop that calculates the chi2
// Main loop that calculates the chi2 - run in parallel using OpenMP
#pragma omp parallel for
			for (long int e = 0; e < emax; e++) {
				double mdass = std::get<0>(vect_2D[e]);
				double mcorr = std::get<1>(vect_2D[e]);
				double like_event = 0.0;
				for (int i = 0; i < Config::ncontr; i++) {
					double md_like, mb_like;
					md_like = D_PDFs_get[i]->EvalPDF(&mdass, &param[Config::nvar_offset_md[i]]);
					mb_like = B_PDFs_get[i]->EvalPDF(
						&mcorr, &param[Config::nvar_all_md + Config::nvar_offset_mb[i]]);
					if (mb_like < 0.0 || mb_like > 1.0 || md_like < 0.0 || md_like > 1.0 ||
						frac[i] < 0.0 || frac[i] > 1.0) {
						chi2 += 1e15;
						continue;
					}
					like_event += md_like * mb_like * frac[i];
				}
				vect_chi2[e] = -TMath::Log(like_event);
				/*
				{
					std::lock_guard<std::mutex> guard(my_mutex);
					std::cout<<"Chi2: "<<likelihood<<std::flush<<std::endl;
				}
				*/
			}
		}

		// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		double chi2_threads = 0.0;
		if (avx) {
			size_t vecsize;
			if (Config::binned)
				vecsize = size_t(hist2D.GetNbinsX());
			else
				vecsize = (size_t)vect_2D.size();
			chi2_threads = fastAccurate<Method::Kahan, 4>(vect_chi2, vecsize);
		} else {
			double sum = 0, c = 0;
			for (long unsigned int i = 0; i < vect_2D.size(); i++) {
				ksum(sum, c, vect_chi2[i]);
			}
			chi2_threads = sum + c;
		}
		delete[] vect_chi2;
		chi2 += chi2_threads;
		// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>
		// (end - begin).count() << "[ns]" << std::endl;
		if (!Config::binned) chi2 += log(double(emax));	 // Extended ML

		return chi2;
	};
	return fchi2;
}

std::function<double(const double*)> wrap_chi2_simultanous(
	const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get,
	const std::vector<std::shared_ptr<PDFInterface>>& B_PDFs_get,
	std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD,
	const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB,
	const std::vector<std::vector<double>>& dMC_MB,
	const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH1D hist_MD, TH1D hist_MB) {
	long int emax = vect_2D.size();
	std::cout << emax << " yield of the sample \n";
	auto fchi2 = [D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB, emax,
				  replaceIndexVect, int_choose_fit, hist_MD, hist_MB](const double* par) -> double {
		// Caclulate chi2
		// Get the pointer in parameters array that correspond to fraction defining parameters????
		const double* pa = &par[Config::nvar_all_md + Config::nvar_all_mb];
		// Calculate fractions
		// The parameters pa[] define the fractions, we have 6 fractions but 5 independent
		// parameters. The parametrisation is arbitrary
		double frac[Config::ncontr];

		double chi2 = 0.0;
		double sum_frac = 0.0;
		for (int i = 0; i < Config::ncontr - 1; i++) {
			frac[i] = abs(pa[i]);
			sum_frac += frac[i];
		}
		frac[Config::ncontr - 1] = abs(1.0 - sum_frac);
		sum_frac += frac[Config::ncontr - 1];
		double tmp = 1.0e10 * (sum_frac - 1.0) * (sum_frac - 1.0);
		if (!Config::binned) tmp *= emax;
		chi2 += tmp;

		// Extract the parameters and add some constraints
		double param[Config::nvar_all_md + Config::nvar_all_mb];
		for (int ivar = 0; ivar < (Config::nvar_all_md + Config::nvar_all_mb); ivar++) {
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var : replaceIndexVect) {
			param[irep_var.first] = par[irep_var.second];
		}

		// Calculate normalisation integrals
		for (int i = 0; i < Config::ncontr; i++) {
			D_PDFs_get[i]->CalcIntegral(&param[Config::nvar_offset_md[i]], Config::minDM,
										Config::maxDM);
			B_PDFs_get[i]->CalcIntegral(&param[Config::nvar_all_md + Config::nvar_offset_mb[i]],
										Config::minBMcorr, Config::maxBMcorr);
		}

		int nbins_md = hist_MD.GetNbinsX();
		int nbins_mb = hist_MB.GetNbinsY();
		if (int_choose_fit !=
			dictionaryChooseFit.at("frac")) {  // don't calculate if only fractions fit
			for (int i = 0; i < Config::ncontr; i++) {
				// TODO set dxx for D comb background to 0 and remove this
				for (int ivar = 0; ivar < Config::nvar_md[i]; ivar++) {
					// Skip the constrain for chebyshev background: {"Chebyshev", 1}
					// if (Config::intshapesDM[i] == 1 || ivar == 1) continue;
					if (dMC_MD[i][ivar] != 0) {
						double tmp = (MC_MD[i][ivar] - param[Config::nvar_offset_md[i] + ivar]) *
									 (MC_MD[i][ivar] - param[Config::nvar_offset_md[i] + ivar]) /
									 (2.0 * (dMC_MD[i][ivar] *
											 dMC_MD[i][ivar]));	 // use the results of MC fits
						double scale = 1.0;
						chi2 += scale * tmp;
					}
				}
			}
		}
		if (int_choose_fit !=
			dictionaryChooseFit.at("frac")) {  // don't calculate if only fractions fit
			for (int i = 0; i < Config::ncontr; i++) {
				for (int ivar = 0; ivar < Config::nvar_mb[i]; ivar++) {
					if (dMC_MB[i][ivar] != 0) {
						double tmp =
							(MC_MB[i][ivar] -
							 param[Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar]) *
							(MC_MB[i][ivar] -
							 param[Config::nvar_all_md + Config::nvar_offset_mb[i] + ivar]) /
							(2.0 *
							 (dMC_MB[i][ivar] * dMC_MB[i][ivar]));	// use the results of MC fits

						double scale = 1.0;
						// if(int_choose_fit == dictionaryChooseFit.at("BM"))
						//	scale = 10000.0;
						chi2 += scale * tmp;
					}
				}
			}
		}

		double* vect_chi2_md =
			new double[vect_2D.size()];	 // Per event results - required to efficiently calculate a
										 // Kahan compensated sum
		double* vect_chi2_mb =
			new double[vect_2D.size()];	 // Per event results - required to efficiently calculate a
										 // Kahan compensated sum
		for (int jj = 0; jj < int(vect_2D.size()); jj++) vect_chi2_md[jj] = 0.0;
		for (int jj = 0; jj < int(vect_2D.size()); jj++) vect_chi2_mb[jj] = 0.0;
		if (Config::binned) {
			double bin_width_mb = (Config::maxBMcorr - Config::minBMcorr) / double(nbins_mb);
			double bin_width_md = (Config::maxDM - Config::minDM) / double(nbins_md);
			double histev_md = hist_MD.Integral();
			double histev_mb = hist_MB.Integral();
#pragma omp parallel for
			for (int bin_md = 1; bin_md <= nbins_md; bin_md++) {
				double sum_contr = 0.0;
				for (int i = 0; i < Config::ncontr; i++) {
					double mdass = hist_MD.GetXaxis()->GetBinCenter(bin_md);
					double md_val =
						D_PDFs_get[i]->EvalPDF(&mdass, &param[Config::nvar_offset_md[i]]);
					sum_contr +=
						histev_md * bin_width_md * frac[i] * md_val;
					// std::cout << sum_contr << "  sumcontr  " << frac[i] << "  " << md_val <<
					// "  " << mb_val << std::endl;
					if (frac[i] < 0.0 || frac[i] > 1.0) {
						chi2 += 1e15;
						continue;
					}
				}
				double cont = (double)hist_MD.GetBinContent(bin_md);
				double err = (double)hist_MD.GetBinError(bin_md);
				if (err != 0.0) {
					vect_chi2_md[bin_md - 1] +=
						0.5 * (sum_contr - cont) * (sum_contr - cont) / err / err;
				}
			}
#pragma omp parallel for
			for (int bin_mb = 1; bin_mb <= nbins_mb; bin_mb++) {
				double sum_contr = 0.0;
				for (int i = 0; i < Config::ncontr; i++) {
					double mcorr = hist_MB.GetYaxis()->GetBinCenter(bin_mb);
					double mb_val = B_PDFs_get[i]->EvalPDF(
						&mcorr, &param[Config::nvar_all_md + Config::nvar_offset_mb[i]]);
					sum_contr +=
						histev_mb * bin_width_mb * frac[i] * mb_val;
					// std::cout << sum_contr << "  sumcontr  " << frac[i] << "  " << md_val <<
					// "  " << mb_val << std::endl;
					if (frac[i] < 0.0 || frac[i] > 1.0) {
						chi2 += 1e15;
						continue;
					}
				}
				double cont = (double)hist_MB.GetBinContent(bin_mb);
				double err = (double)hist_MB.GetBinError(bin_mb);
				if (err != 0.0) {
					vect_chi2_mb[bin_mb - 1] +=
						0.5 * (sum_contr - cont) * (sum_contr - cont) / err / err;
				}
			}
		} else {

			// Main loop that calculates the chi2
// Main loop that calculates the chi2 - run in parallel using OpenMP
#pragma omp parallel for
			for (long int e = 0; e < emax; e++) {
				double mdass = std::get<0>(vect_2D[e]);
				double mcorr = std::get<1>(vect_2D[e]);
				double like_event_md = 0.0;
				double like_event_mb = 0.0;
				for (int i = 0; i < Config::ncontr; i++) {
					double md_like, mb_like;
					md_like = D_PDFs_get[i]->EvalPDF(&mdass, &param[Config::nvar_offset_md[i]]);
					mb_like = B_PDFs_get[i]->EvalPDF(
						&mcorr, &param[Config::nvar_all_md + Config::nvar_offset_mb[i]]);
					if (mb_like < 0.0 || mb_like > 1.0 || md_like < 0.0 || md_like > 1.0 ||
						frac[i] < 0.0 || frac[i] > 1.0) {
						chi2 += 1e15;
						continue;
					}
					like_event_md += md_like * frac[i];
					like_event_mb += mb_like * frac[i];
				}
				vect_chi2_md[e] = -TMath::Log(like_event_md);
				vect_chi2_mb[e] = -TMath::Log(like_event_mb);
				/*
				{
					std::lock_guard<std::mutex> guard(my_mutex);
					std::cout<<"Chi2: "<<likelihood<<std::flush<<std::endl;
				}
				*/
			}
		}

		// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		double chi2_threads = 0.0;
		if (avx) {
			size_t vecsize_md;
			size_t vecsize_mb;
			if (Config::binned) {
				vecsize_md = size_t(hist_MD.GetNbinsX());
				vecsize_mb = size_t(hist_MB.GetNbinsX());
			} else {
				vecsize_md = (size_t)vect_2D.size();
				vecsize_mb = (size_t)vect_2D.size();
			}
			chi2_threads = fastAccurate<Method::Kahan, 4>(vect_chi2_md, vecsize_md);
			chi2_threads += fastAccurate<Method::Kahan, 4>(vect_chi2_mb, vecsize_mb);
		} else {
			double sum = 0, c = 0;
			for (long unsigned int i = 0; i < vect_2D.size(); i++) {
				ksum(sum, c, vect_chi2_md[i]);
			}
			chi2_threads = sum + c;
			sum = c = 0;
			for (long unsigned int i = 0; i < vect_2D.size(); i++) {
				ksum(sum, c, vect_chi2_mb[i]);
			}
			chi2_threads += sum + c;
		}
		delete[] vect_chi2_md;
		delete[] vect_chi2_mb;
		chi2 += chi2_threads;
		// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds>
		// (end - begin).count() << "[ns]" << std::endl;

		// TODO is this now correct ?? Should it be 2 x emax ??
		if (!Config::binned) chi2 += log(double(emax));	 // Extended ML

		return chi2;
	};
	return fchi2;
}

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
