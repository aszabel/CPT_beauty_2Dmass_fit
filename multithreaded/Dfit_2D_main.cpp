#include <string>
#include <fstream>
#include <iomanip>
#include "omp.h"
#include "json.hpp"
#include "config.h"

// ROOT includes
#include "TROOT.h"
#include "TH2D.h"
#include "TChain.h"
#include "TString.h"
#include "TF2.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/ProbFuncMathCore.h"

// CPT_beauty_2Dmass_fit
#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "FastSum.h"
#include "ChebyshevPDF.h"

#include <functional>
#include <memory>

using json = nlohmann::json;
using namespace cpt_b0_analysis;

#include <mutex>
std::mutex my_mutex;

// bool avx = false;
bool avx = true;

/**
 * @brief Draw the results of a 2D fit
 *
 * @param res
 * @param nevents
 * @param name
 * @param sign
 */
// At the moment this is not used
// void Draw(const double *res, double nevents, TString name[], int sign);

/**
 * @brief A 2D fitter for ..
 *
 * @param argc
 * @param argv
 * @return int
 */
std::function<double(const double*)> wrap_chi2(const std::vector<PDFInterface*>& D_PDFs_get, const std::vector<PDFInterface*>& B_PDFs_get, std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB, const std::vector<std::vector<double>>& dMC_MB);



int main(int argc, char *argv[])
{
	// The first, fundamental operation to be performed in order to make ROOT
	// thread-aware.
	ROOT::EnableThreadSafety();

	// Increase printout precision
	// std::cout<<std::setprecision(40);

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
        int sign = Config::sign;
        std::string input_file = Config::input_file;

        int nentries = Config::nentries;

        double tolerance = Config::tolerance;

        double muPTmin = Config::muPTmin;
        double muPmin = Config::muPmin;
        static double eta_min = Config::eta_min;
        static double eta_max = Config::eta_max;

        //Global
        double minDM = Config::minDM;
        double maxDM = Config::maxDM;
        double minBMcorr = Config::minBMcorr;
        double maxBMcorr = Config::maxBMcorr;
        const int nvar_md = 7;//Config::nvar_md;
        const int nvar_mb = 7;///usr/bin/../lib/gcc/x86_64-linux-gnu/11/../../../../include/c++/11/bits/stl_uninitialized.hConfig::nvar_mb;
        const int ncontr = 6;//Config::ncontr;
        const int nbins = Config::nbins;
	std::vector<std::string> fixVect = Config::fixVect;


	const std::vector<std::vector<double>>& MC_MD = Config::MC_MD;
	const std::vector<std::vector<double>>& dMC_MD = Config::dMC_MD;
	const std::vector<std::vector<double>>& MC_MB = Config::MC_MB;
	const std::vector<std::vector<double>>& dMC_MB = Config::dMC_MB;

	
	std::string minName = "Minuit2";
	std::string algoName = "";
	ROOT::Math::Minimizer *min =
		ROOT::Math::Factory::CreateMinimizer(minName, algoName);

	// Set tolerance , etc...
	// TODO load those values from a config file
	min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
	min->SetMaxIterations(10000);		// for GSL
	min->SetTolerance(tolerance);
	min->SetPrintLevel(2);
	// min->SetStrategy(2);
	// min->SetPrecision(0.00001);

	// Load data set
	TChain ch("BlindedTree");
	// TODO load the data filename from config file
	ch.Add(input_file.c_str());
	double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
	bool charge;

	std::vector<std::pair<double, double>> vect_2D;
	TH2D *hist = new TH2D("hist", "", nbins, minDM, maxDM, nbins, minBMcorr, maxBMcorr);
	ch.SetBranchAddress("B_M", &B_M);
	ch.SetBranchAddress("missPT", &missPT);
	ch.SetBranchAddress("D_M", &D_M);
	ch.SetBranchAddress("mu_PT", &mu_PT);
	ch.SetBranchAddress("mu_P", &mu_P);
	ch.SetBranchAddress("mu_eta", &mu_eta);
	ch.SetBranchAddress("K_PT", &K_PT);
	ch.SetBranchAddress("truecharge", &charge);

	if (nentries < 0) nentries = ch.GetEntries(); 
	if (nentries > ch.GetEntries()){
		std::cerr << "The value of 'nentries' exceeds the number of events in the file." << std::endl;
                return 1;
        }
	std::cout << nentries << " number of entries read"  <<std::endl;
	for (int i=0; i<nentries; ++i){
		ch.GetEntry(i);
		if (mu_PT < muPTmin || mu_P < muPmin || mu_eta < eta_min || mu_eta > eta_max)
			continue;
		double B_MMcorr = B_M + 2.0 * missPT;
		if (B_MMcorr < minBMcorr || B_MMcorr > maxBMcorr)
			continue;
		if (D_M < minDM || D_M > maxDM)
			continue;
		if (int(charge) == sign)
		{
			vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
			hist->Fill(D_M, B_MMcorr);
		}
	}


	// Initial fit parameter values taken from 1D fits to MC and Side Bands

		

	// Define Minuit fit variables for M_D

	const auto& contrName = Config::contrName;
	const auto& varname_md = Config::varname_md;

	double starting_point[(nvar_md+nvar_mb)*ncontr+ncontr-1];
	for (int i = 0; i < ncontr; i++)
	{
		for (int ivar = 0; ivar < nvar_md; ivar++)
		{
			min->SetVariable(i * nvar_md + ivar, (TString::Format("%s_%s", contrName[i].c_str(), varname_md[ivar].c_str())).Data(), MC_MD[i][ivar], dMC_MD[i][ivar] + 1.0e-11);
			starting_point[i * nvar_md + ivar] = MC_MD[i][ivar];
			// min->FixVariable(i*nvar_md+ivar);
			//min->SetVariableLowerLimit(i * nvar_md + ivar, 1e-13);
			// TODO load limits and fix values from a config file
			if (ivar == 3)
				min->SetVariableLimits(i * nvar_md + ivar, -1.0, 1.0);

			if (ivar == 5 || ivar == 6)
				min->FixVariable(i * nvar_md + ivar);
			if (i == 4 && ivar != 0 && ivar != 1)
				min->FixVariable(i * nvar_md + ivar); // combinatorial
			// if (i==4) min->FixVariable(i*nvar_md+ivar); // combinatorial
			// if(ivar!=1&&ivar!=0)min->FixVariable(i*nvar_md+ivar);
			if (i != 0 && ivar == 1)
				min->FixVariable(i * nvar_md + ivar);
		}
	}

	// Define Minuit fit variables for M_B
	const auto & varname_mb = Config::varname_mb;

	for (int i = 0; i < ncontr; i++)
	{
		// TODO load limits and fixed value from the config file
		for (int ivar = 0; ivar < nvar_mb; ivar++)
		{
			min->SetVariable(ncontr * nvar_md + i * nvar_mb + ivar, (TString::Format("%s_%s", contrName[i].c_str(), varname_mb[ivar].c_str())).Data(), MC_MB[i][ivar], dMC_MB[i][ivar] + 1.0e-11);
			starting_point[ncontr * nvar_md + i * nvar_mb + ivar] = MC_MB[i][ivar];
			//min->FixVariable(ncontr*nvar_md+i*nvar_mb+ivar);
			//min->SetVariableLowerLimit(ncontr * nvar_md + i * nvar_mb + ivar, 1e-13);
			if (ivar == 2)
				min->SetVariableLimits(ncontr * nvar_md + i * nvar_mb + ivar, -1.0, 1.0);
		}
	}

	// Set initial fraction values
	const auto&fracInit = Config::fracInit;

	for (int i = 0; i < ncontr - 1; i++)
	{
		min->SetVariable((nvar_md + nvar_mb) * ncontr + i, (TString::Format("par_frac%d", i)).Data(), fracInit[i], 0.01);
		starting_point[(nvar_md + nvar_mb) * ncontr + i] = fracInit[i];
	}


	// Q: Define ??????
	// Define the error setimation parameter in minuit for 1 sigma and ncontr -1 free parameters
	double CL_normal = ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1); //1 sigma ~68%
	min->SetErrorDef(TMath::ChisquareQuantile(CL_normal, ncontr-1));// ncontr-1 free fraction parameters, other parameters have gaussian contraints base on MC fits.

	const auto& D_PDFs = Config::getVectorPDFs("Dmass");
	if (D_PDFs.size()!=ncontr){
		std::cout<< " NO D_PDFs \n";
		return 1;//std::function<double(const double*)>{};
	}
	for (int i=0; i<ncontr; ++i){
		auto pdf = D_PDFs[i].get();
		std::cout << pdf << std::endl;

		if (!pdf){
			std::cout<< "Nullptr passed as pdf\n";
			return 1;//std::function<double(const double*)>{};
		}
	}
		
	const auto& B_PDFs = Config::getVectorPDFs("Bmass");
        if (B_PDFs.size()!=ncontr){
                std::cout<< " NO B_PDFs \n";    
                return 1;//std::function<double(const double*)>{};
        }
        for (int i=0; i<ncontr; ++i){
                auto pdf = B_PDFs[i].get();
		std::cout << pdf << std::endl;
                if (!pdf){
                        std::cout<< "Nullptr passed as pdf\n";
                        return 1;//std::function<double(const double*)>{};
                }
        }
	std::vector<PDFInterface *> D_PDFs_get;
	std::vector<PDFInterface *> B_PDFs_get;
	for (int i=0; i<ncontr; ++i){
		D_PDFs_get.push_back(D_PDFs[i].get());
		B_PDFs_get.push_back(B_PDFs[i].get());
	}
	//staring point for stability checks
	for (int ivar=0; ivar<(nvar_md+nvar_mb)*ncontr+ncontr-1; ivar++)
		min->SetVariableValue(ivar, starting_point[ivar]);

	//list of fixed variables form config
	for (const auto& fix: fixVect){
		min->FixVariable(min->VariableIndex(fix));
		//std::cout << fix << "  " << min->VariableIndex(fix) << std::endl;
	}

	// Start the minimization
	// Define a fit function for Minuit
	auto fchi2 = wrap_chi2(D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB);
	ROOT::Math::Functor f(fchi2, (nvar_md + nvar_mb) * ncontr + ncontr - 1);
	min->SetFunction(f);
	min->Minimize();

	// Print the fit results
	std::ofstream results(TString::Format("results_%d.txt", sign));
	for (int i = 0; i < (nvar_md + nvar_mb) * ncontr + ncontr - 1; i++)
	{
		results << min->X()[i] << "  " << min->Errors()[i] << std::endl;
	}
	results.close();
	return 0;
}
std::function<double(const double*)> wrap_chi2(const std::vector<PDFInterface*>& D_PDFs_get, const std::vector<PDFInterface*>& B_PDFs_get, std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB, const std::vector<std::vector<double>>& dMC_MB){

        
	
        const int ncontr = Config::ncontr;

	
	


	double minDM = Config::minDM;
        double maxDM = Config::maxDM;
        double minBMcorr = Config::minBMcorr;
        double maxBMcorr = Config::maxBMcorr;
        const int nvar_md = Config::nvar_md;
        const int nvar_mb = Config::nvar_mb;




	long int emax = vect_2D.size();
	std::cout << emax << "emax \n";
	auto fchi2 = [D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB, emax, ncontr, minDM, maxDM, minBMcorr, maxBMcorr, nvar_md, nvar_mb](const double *par) -> double 
	{

	// Caclulate chi2
			// Get the pointer in parameters array that correspond to fraction defining parameters????
		const double *pa = &par[ncontr * (nvar_md + nvar_mb)];
		// Calculate fractions
		// Q: Why do we recalculate fractions? What does the minuti minimize - what is stored in `pa` ???
		// The parameters pa[] define the fractions, we have 6 fractions but 5 independent parameters. The parametrisation is arbitrary
		double frac[ncontr];
		frac[0] = 1 - abs(pa[0]);
		frac[1] = abs(pa[0]) * (1.0 - abs(pa[1]));
		frac[2] = abs(pa[0]) * abs(pa[1]) * (1.0 - abs(pa[2]));
		frac[3] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * (1.0 - abs(pa[3]));
		frac[4] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * abs(pa[3]) * (1.0 - abs(pa[4]));
		frac[5] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * abs(pa[3]) * abs(pa[4]);
		// Extract the parameters and add some constraints
		double param[ncontr * (nvar_md + nvar_mb)];
		for (int i = 0; i < ncontr; i++)
		{
			for (int ivar = 0; ivar < nvar_md; ivar++)
			{
				param[i * nvar_md + ivar] = par[i * nvar_md + ivar];
				// TODO find a better way ...
				// Define this in a json config
				if (ivar == 1 && i != 2 && i != 4)
					param[i * nvar_md + ivar] = par[1];
			}
			for (int ivar = 0; ivar < nvar_mb; ivar++)
			{
				int ipar = ncontr * nvar_md + i * nvar_mb + ivar;
				param[ipar] = par[ipar];
			}
		}
		/*
		for (int i=0; i<ncontr; ++i){
			std::cout<<D_PDFs[i].get() << "  " << B_PDFs[i].get() << std::endl;
		}
		*/
		// Calculate normalisation integrals
		for (int i = 0; i < ncontr; i++)
		{
			D_PDFs_get[i]->CalcIntegral(&param[i * nvar_md], minDM, maxDM);
			B_PDFs_get[i]->CalcIntegral(&param[ncontr * nvar_md + i * nvar_mb], minBMcorr, maxBMcorr);
		}

		// Main loop that calculates the chi2
		double chi2 = 0.0;
		double *vect_chi2 = new double[vect_2D.size()]; // Per event results - required to efficiently calculate a Kahan compensated sum

// Main loop that calculates the chi2 - run in parallel using OpenMP
#pragma omp parallel for
		for (long int e = 0; e < emax; e++)
		{
			double mdass = std::get<0>(vect_2D[e]);
			double mcorr = std::get<1>(vect_2D[e]);
			double likelihood = 0.0;
			for (int i = 0; i < ncontr; i++)
			{
				double md_like, mb_like;
				md_like = D_PDFs_get[i]->EvalPDF(&mdass, &param[i * nvar_md]);
				mb_like = B_PDFs_get[i]->EvalPDF(&mcorr, &param[ncontr * nvar_md + i * nvar_mb]);
				likelihood += md_like * mb_like * frac[i];
			}
			if (likelihood > 0.0)
				vect_chi2[e] = -2.0 * log(likelihood);
			/*
			{
				std::lock_guard<std::mutex> guard(my_mutex);
				std::cout<<"Chi2: "<<likelihood<<std::flush<<std::endl;
			}
			*/
		}

		// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		if (avx)
		{
			chi2 = fastAccurate<Method::Kahan, 4>(vect_chi2, (size_t)vect_2D.size());
		}
		else
		{
			double sum = 0, c = 0;
			for (long unsigned int i = 0; i < vect_2D.size(); i++)
			{
				ksum(sum, c, vect_chi2[i]);
			}
			chi2 = sum + c;
		}
		delete [] vect_chi2;
		// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;

		// Add additional constraints based on the templates from 1D fits
		for (int i = 0; i < ncontr; i++)
		{
			// TODO set dxx for D comb background to 0 and remove this
			for (int ivar = 0; ivar < nvar_md; ivar++)
			{
				if (dMC_MD[i][ivar] != 0)
					chi2 += (MC_MD[i][ivar] - param[i * nvar_md + ivar]) * (MC_MD[i][ivar] - param[i * nvar_md + ivar]) / (dMC_MD[i][ivar] * dMC_MD[i][ivar]); // use the results of MC fits
			}
		}
		for (int i = 0; i < ncontr; i++)
		{
			for (int ivar = 0; ivar < nvar_mb; ivar++)
			{
				if (dMC_MB[i][ivar] != 0)
					chi2 += (MC_MB[i][ivar] - param[ncontr * nvar_md + i * nvar_mb + ivar]) * (MC_MB[i][ivar] - param[ncontr * nvar_md + i * nvar_mb + ivar]) / (dMC_MB[i][ivar] * dMC_MB[i][ivar]); // use the results of MC fits
			}
		}

		// std::cerr<<"Chi2: "<<chi2<<std::endl;
		return chi2;
	};
	return fchi2;
}

