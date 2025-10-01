/** TODO TODO
 * - Remove all contributions besides sidebands - they are not used anyway
 * - Set output directory
 * - Add chebyshev of higher orders ...
 * - Test unbinned version
 */

#include <fstream>
#include <iomanip>
#include <string>

#include "config.h"
#include "json.hpp"
#include "omp.h"

// ROOT includes
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/ProbFuncMathCore.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TMath.h"
#include "TPad.h"
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
#include "TH1D.h"

// using json = nlohmann::json;
using namespace cpt_b0_analysis;

#include <mutex>
std::mutex my_mutex;

// bool avx = false;
bool avx = true;

/**
 * @brief Wrapper function for chi2 calculation
 *
 * @param D_PDFs_get vector of D_M PDF instances
 * @param data D_M and B_M mass values for considered events
 * @param MC_MD
 * @param dMC_MD
 * @param replaceIndexVect set of rules for replacement of fit parameter values,
 *                         used to obtain common values between several PDF components
 * @param int_choose_fit fit type: "frac", "BM", "DM+BMfixed", "all"
 * @param hist1D
 * @return int
 */
std::function<double(const double*)> wrap_chi2(
	const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get, std::vector<double> data,
	const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD,
	const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH1D hist1D);

/**
 * @brief Draw pull histogram
 *
 * @param c Canvas to draw on
 * @param hist data
 * @param tf1_sum PDF to compare with data
 * @return int
 */
void Draw_pull(TCanvas* c, TH1D hist, TF1* tf1_sum);

/**
 * @brief Fit combinatorial fraction from D_M sidebands by fixing all other parameters
 *
 * @param argc
 * @param argv - takes a single argument - path to a config file
 * @return int
 */
int main(int argc, char* argv[]) {
	// The first, fundamental operation to be performed in order to make ROOT
	// thread-aware.
	// ROOT::EnableThreadSafety();

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
	int nentries = Config::nentries;

	// Load data set - only single data set is used
	TChain ch((Config::chainName).c_str());
	ch.Add(Config::input_files[0].c_str());
	double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
	double B_MMcorr;
	bool charge;
	int frac_index = 100;
	const int nbins = 100;

	// D_M histogram for sidebands
	TH1D hist1D("hist1D", "", nbins, Config::minDM, Config::maxDM);
	double bin_width = (Config::maxDM - Config::minDM) / double(nbins);

	std::vector<double> data = {};
	if (!Config::isMC) {
		ch.SetBranchAddress("B_M", &B_M);
		ch.SetBranchAddress("missPT", &missPT);
		ch.SetBranchAddress("D_M", &D_M);
		ch.SetBranchAddress("mu_PT", &mu_PT);
		ch.SetBranchAddress("mu_P", &mu_P);
		ch.SetBranchAddress("mu_eta", &mu_eta);
		ch.SetBranchAddress("K_PT", &K_PT);
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
		if (!Config::isMC) {
			if (mu_PT < Config::muPTmin || mu_P < Config::muPmin || mu_eta < Config::eta_min ||
				mu_eta > Config::eta_max)
				continue;
			if (int(charge) != Config::sign) continue;
			B_MMcorr = B_M + 2.0 * missPT;
		}
		if (B_MMcorr < Config::minBMcorr || B_MMcorr > Config::maxBMcorr) continue;

		if (D_M < Config::minDM || D_M > Config::maxDM) continue;
		data.push_back(D_M);

		// Histogram only D_M sidebands
		if (D_M < hist1D.GetBinLowEdge(hist1D.FindBin(1830)) ||
			D_M >= hist1D.GetBinLowEdge(hist1D.FindBin(1910)))
			hist1D.Fill(D_M);

		// TODO_DOCS Calc number of entries ????
		if (Config::isMC) frac_indeces[frac_index]++;
	}

	// Initial fit parameter values are taken from 1D fits to MC and Side Bands

	// Define Minuit fit variables for M_D
	const int n_all = Config::nvar_all_md + Config::ncontr;
	// Get index of the sidebands contribution
	const int sb_idx =
		std::distance(Config::contrName.begin(),
					  std::find(Config::contrName.begin(), Config::contrName.end(), "sidebands"));
	bool previous_fit = false;
	bool start_scratch = true;
	double starting_point[n_all];
	std::string minName = "Minuit2";
	std::string algoName = "";
	int loop_count = 0;
	for (auto& int_choose_fit : Config::int_choose_fits) {
		ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

		// Set tolerance , etc...
		min->SetMaxFunctionCalls(Config::functionCalls);  // for Minuit/Minuit2
		min->SetTolerance(Config::tolerance[loop_count]);
		loop_count++;
		min->SetPrintLevel(Config::printLevel);
		// min->SetStrategy(2);

		// Define Minuit fit variables for M_D
		for (int i = 0; i < Config::ncontr; i++) {
			for (int ivar = 0; ivar < Config::nvar_md[i]; ivar++) {
				// Set initial variable values from 1D fit
				// Set step size based on 1D fit uncertainty
				min->SetVariable(Config::nvar_offset_md[i] + ivar,
								 (TString::Format("%s_%s", Config::contrName[i].c_str(),
												  Config::varname_md[i][ivar].c_str()))
									 .Data(),
								 Config::MC_MD[i][ivar], Config::dMC_MD[i][ivar] + 1.0e-11);
				if (i != sb_idx) {
					// Fix parameters for all contributions besides sidebands
					min->FixVariable(Config::nvar_offset_md[i] + ivar);
				} else {
					if (ivar >= Config::n_sideband)
						min->FixVariable(Config::nvar_offset_md[i] + ivar);
				}
				// Read initial value from 1D fit for "from scratch" fit
				if (start_scratch)
					starting_point[Config::nvar_offset_md[i] + ivar] = Config::MC_MD[i][ivar];
			}
		}

		// Set Limits on variables from config file
		for (auto it = Config::varLimitsMap.begin(); it != Config::varLimitsMap.end(); ++it) {
			int index_var = min->VariableIndex(it->first);
			if (index_var == -1) {
				std::cerr << "Error in limiting parameters: param " << it->first << " not found.\n";
				return 1;
			}
			auto pairlims = it->second;
			min->SetVariableLimits(index_var, pairlims.first, pairlims.second);
		}

		for (int i = 0; i < Config::ncontr; i++) {
			// TODO if this is 1D M_D fit why are we setting fractions after M_B params ???
			min->SetVariable(Config::nvar_all_md + i, (TString::Format("par_frac%d", i)).Data(),
							 Config::fracInit[i], 0.001);
			if (i != sb_idx) {
				min->SetVariableValue(Config::nvar_all_md + i, 0.0);
				min->FixVariable(Config::nvar_all_md + i);
			}
			if (start_scratch) starting_point[Config::nvar_all_md + i] = Config::fracInit[i];
		}

		// Define the error setimation parameter in minuit for 1 sigma and ncontr -1 free parameters
		double CL_normal = ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1);	// 1 sigma ~68%
		min->SetErrorDef(TMath::ChisquareQuantile(
			CL_normal, Config::ncontr - 1));  // ncontr-1 free fraction parameters, other parameters
											  // have gaussian contraints base on MC fits.

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

		// list of fixed variables form config*/
		for (const auto& fix : Config::fixVect) {
			min->FixVariable(min->VariableIndex(fix));
			std::cout << "Fix: " << fix << "  " << min->VariableIndex(fix) << std::endl;
		}

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
		auto fchi2 = wrap_chi2(D_PDFs, data, Config::MC_MD, Config::dMC_MD, replaceIndexVect,
							   int_choose_fit, hist1D);
		ROOT::Math::Functor f(fchi2, n_all);
		min->SetFunction(f);
		min->Minimize();

		// Print the fit results
		std::ofstream results(
			TString::Format("results1D_DM_%d_%d.txt", Config::sign, int_choose_fit));
		results << min->Status() << "  " << min->MinValue() << std::endl;
		if (min->Status() != 0) {
			std::cout << "Bad status of fit " << min->Status() << std::endl;
			previous_fit = false;
		} else {
			previous_fit = true;
			start_scratch = false;
		}
		const double* pa = &min->X()[Config::nvar_all_md];
		for (int i = 0; i < Config::nvar_mb[sb_idx]; i++) {
			results << min->X()[Config::nvar_offset_md[sb_idx] + i] << "  "
					<< min->Errors()[Config::nvar_offset_md[sb_idx] + i] << std::endl;
		}
		results << min->X()[Config::nvar_all_md + sb_idx] << "  "
				<< min->Errors()[Config::nvar_all_md + sb_idx] << std::endl;

		results.close();

		for (int i = 0; i < Config::ncontr; i++)
			std::cout << abs(pa[i]) << "  frac" << i << "  "
					  << double(frac_indeces[i]) / double(data.size()) << std::endl;

		D_PDFs[sb_idx].get()->CalcIntegral(&min->X()[Config::nvar_offset_md[sb_idx]],
										   hist1D.GetBinLowEdge(hist1D.FindBin(1830)),
										   hist1D.GetBinLowEdge(hist1D.FindBin(1910)));
		double integral_signal_range = D_PDFs[sb_idx].get()->getIntegral();
		std::cout << integral_signal_range << " intsig\n";
		D_PDFs[sb_idx].get()->CalcIntegral(&min->X()[Config::nvar_offset_md[sb_idx]], Config::minDM,
										   Config::maxDM);
		double integral_full_range = D_PDFs[sb_idx].get()->getIntegral();
		std::cout << integral_full_range << " intfull\n";

		double nevents = (double)hist1D.Integral();
		TF1* sidebandtf1 = new TF1(
			"stf1",
			[D_PDFs, nevents, integral_full_range, integral_signal_range, bin_width, sb_idx](
				double* x, double* par) -> double {
				return bin_width * nevents * integral_full_range /
					   (integral_full_range - integral_signal_range) *
					   D_PDFs[sb_idx].get()->EvalPDF(x, par);
			},
			Config::minDM, Config::maxDM, Config::n_sideband);
		sidebandtf1->SetParameters(&min->X()[Config::nvar_offset_md[sb_idx]]);

		std::cout << sidebandtf1->Integral(Config::minDM, Config::maxDM) << " intside "
				  << sidebandtf1->Integral(hist1D.GetBinLowEdge(hist1D.FindBin(1830)),
										   hist1D.GetBinLowEdge(hist1D.FindBin(1910)))
				  << " nevents " << nevents << std::endl;

		TCanvas* c = new TCanvas("x", "", 500, 500);
		Draw_pull(c, hist1D, sidebandtf1);
		c->SaveAs("sideband.pdf");
		if (min) {
			delete min;
		}
	}
	return 0;
}

void Draw_pull(TCanvas* c, TH1D hist, TF1* tf1_sum) {
	c->cd();
	TPad* pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
	pad1->SetLogy();
	pad1->Draw();
	pad1->cd();

	hist.SetStats(kFALSE);
	hist.SetMinimum(1.0);
	hist.DrawClone("ep");
	tf1_sum->SetLineColor(kBlack);
	tf1_sum->DrawClone("same");

	// hist->Sumw2();
	TH1D histpull1D(hist);
	histpull1D.SetMinimum();
	for (int bin = 1; bin <= hist.GetNbinsX(); bin++) {
		double err = hist.GetBinError(bin);
		if (err == 0.0) continue;
		double diff = hist.GetBinContent(bin) - tf1_sum->Eval(hist.GetBinCenter(bin));
		histpull1D.SetBinContent(bin, diff / err);
	}
	c->cd();
	TPad* pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
	pad2->Draw();
	pad2->cd();
	histpull1D.SetStats(kFALSE);
	histpull1D.SetFillColor(kBlue);
	histpull1D.GetYaxis()->SetLabelSize(0.1);
	histpull1D.GetXaxis()->SetLabelSize(0.1);
	histpull1D.DrawClone("hist");
}

std::function<double(const double*)> wrap_chi2(
	const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get, std::vector<double> data,
	const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD,
	const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH1D hist1D) {
	long int emax = data.size();
	std::cout << emax << " yield of the sample \n";
	auto fchi2 = [D_PDFs_get, data, MC_MD, dMC_MD, emax, replaceIndexVect, int_choose_fit,
				  hist1D](const double* par) -> double {
		// Caclulate chi2
		// Get the pointer in parameters array that correspond to fraction defining parameters????
		const double* pa = &par[Config::nvar_all_md];
		// Calculate fractions
		// Q: Why do we recalculate fractions? What does the minuti minimize - what is stored in
		// `pa` ??? The parameters pa[] define the fractions, we have 6 fractions but 5 independent
		// parameters. The parametrisation is arbitrary
		double frac[Config::ncontr];

		double chi2 = 0.0;
		// double sum_frac = 0.0;
		for (int i = 0; i < Config::ncontr; i++) {
			frac[i] = abs(pa[i]);
			// sum_frac+=frac[i];
		}
		// Extract the parameters and add some constraints
		double param[Config::nvar_all_md];
		for (int ivar = 0; ivar < Config::nvar_all_md; ivar++) {
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var : replaceIndexVect) {
			param[irep_var.first] = par[irep_var.second];
		}

		// Calculate normalisation integrals
		for (int i = 0; i < Config::ncontr; i++) {
			D_PDFs_get[i]->CalcIntegral(&param[Config::nvar_offset_md[i]], Config::minDM,
										Config::maxDM);
		}
		double* vect_chi2 = new double[data.size()];  // Per event results - required to efficiently
													  // calculate a Kahan compensated sum
		for (int jj = 0; jj < int(data.size()); jj++) vect_chi2[jj] = 0.0;
		// double test_chi2 = 0.0;
		if (Config::binned) {
			int nbins_md = hist1D.GetNbinsX();
			double bin_width_md = (Config::maxDM - Config::minDM) / double(nbins_md);
			// double histev = hist1D.Integral();
			for (int bin_md = 1; bin_md <= nbins_md; bin_md++) {
				double sum_contr = 0.0;
				for (int i = 0; i < Config::ncontr; i++) {
					double mdass = hist1D.GetXaxis()->GetBinCenter(bin_md);
					double md_val =
						D_PDFs_get[i]->EvalPDF(&mdass, &param[Config::nvar_offset_md[i]]);
					sum_contr += emax * bin_width_md * frac[i] * md_val;
					// std::cout << sum_contr << "  sumcontr  " << frac[i] << "  " << md_val << "  "
					// << mb_val << std::endl;
				}
				double cont = (double)hist1D.GetBinContent(bin_md);
				double err = (double)hist1D.GetBinError(bin_md);
				if (err != 0.0) {
					vect_chi2[bin_md - 1] +=
						0.5 * (sum_contr - cont) * (sum_contr - cont) / err / err;
					// test_chi2 += 0.5 * (sum_contr - cont) * (sum_contr - cont) / err / err;
					//  std::cout << bin_md << "   " << cont << "  " << sum_contr
					//  << "  " << sum_contr/cont << std::endl;
				}
				// std::cout << test_chi2 << " chi2 \n";
				//}
			}
		} else {
			// Main loop that calculates the chi2
			// Main loop that calculates the chi2 - run in parallel using OpenMP
			for (long int e = 0; e < emax; e++) {
				double mdass = data[e];
				double like_event = 0.0;
				for (int i = 0; i < Config::ncontr; i++) {
					double md_like;
					md_like = D_PDFs_get[i]->EvalPDF(&mdass, &param[Config::nvar_offset_md[i]]);
					if (md_like < 0.0 || md_like > 1.0 || frac[i] < 0.0 || frac[i] > 1.0) {
						chi2 += 1e15;
						continue;
					}
					like_event += md_like * frac[i];
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
				vecsize = size_t(hist1D.GetNbinsX());
			else
				vecsize = (size_t)data.size();
			chi2_threads = fastAccurate<Method::Kahan, 4>(vect_chi2, vecsize);
		} else {
			double sum = 0, c = 0;
			for (long unsigned int i = 0; i < data.size(); i++) {
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

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
