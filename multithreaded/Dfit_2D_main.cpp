#include <string>
#include <fstream>
#include <iomanip>
#include "omp.h"
#include "json.hpp"
#include "config.h"

// ROOT includes
#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom.h"
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
#include <tuple>

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
std::function<double(const double*)> wrap_chi2(const std::vector<PDFInterface*>& D_PDFs_get, const std::vector<PDFInterface*>& B_PDFs_get, std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB, const std::vector<std::vector<double>>& dMC_MB, const std::vector<std::pair<int, int>>& replaceIndexVect);



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
        const int nvar_md = Config::nvar_md;
        const int nvar_mb = Config::nvar_mb;
        const int ncontr = Config::ncontr;
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
	min->SetMaxFunctionCalls(Config::functionCalls); // for Minuit/Minuit2
	min->SetTolerance(tolerance);
	min->SetPrintLevel(Config::printLevel);
	// min->SetStrategy(2);

	// Load data set
	TChain ch((Config::chainName).c_str());
	ch.Add(input_file.c_str());
	double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
	bool charge;

	std::vector<std::pair<double, double>> vect_2D;
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

		}
	}

	// Define Minuit fit variables for M_B
	const auto & varname_mb = Config::varname_mb;

	for (int i = 0; i < ncontr; i++)
	{
		for (int ivar = 0; ivar < nvar_mb; ivar++)
		{
			min->SetVariable(ncontr * nvar_md + i * nvar_mb + ivar, (TString::Format("%s_%s", contrName[i].c_str(), varname_mb[ivar].c_str())).Data(), MC_MB[i][ivar], dMC_MB[i][ivar] + 1.0e-11);
			starting_point[ncontr * nvar_md + i * nvar_mb + ivar] = MC_MB[i][ivar];
		}
	}

	//Set Limits on variables
	const std::vector<std::tuple<std::string, double, double>>& varLimitsVect = Config::varLimitsVect;
        for (const auto& varLimits: varLimitsVect){
		int index_var = min->VariableIndex(std::get<0>(varLimits));
                if (index_var == -1 ){
                	std::cerr<< "Error in substituting parameters param " << std::get<0>(varLimits) << " not found.\n";
			return 1;
                }
		min->SetVariableLimits(index_var, std::get<1>(varLimits), std::get<2>(varLimits));
	}


	// Set initial fraction values
	const auto& fracInit = Config::fracInit;
	//auto fracInit = Config::fracInit;
	//double sum = 0.0;
	//for (auto frac: fracInit)
	//	sum += frac;
	//fracInit.push_back(1.0-sum);

	for (int i = 0; i < ncontr - 1; i++)
	//for (int i = 0; i < ncontr; i++)
	{
		min->SetVariable((nvar_md + nvar_mb) * ncontr + i, (TString::Format("par_frac%d", i)).Data(), fracInit[i], 0.01);
		starting_point[(nvar_md + nvar_mb) * ncontr + i] = fracInit[i];
	}


	// Define the error setimation parameter in minuit for 1 sigma and ncontr -1 free parameters
	double CL_normal = ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1); //1 sigma ~68%
	min->SetErrorDef(TMath::ChisquareQuantile(CL_normal, ncontr-1));// ncontr-1 free fraction parameters, other parameters have gaussian contraints base on MC fits.

	const auto& D_PDFs = Config::getVectorPDFs("Dmass");
	if (int(D_PDFs.size())!=ncontr){
		std::cout<< " NO D_PDFs \n";
		return 1;
	}
	for (int i=0; i<ncontr; ++i){
		auto pdf = D_PDFs[i].get();

		if (!pdf){
			std::cout<< "Nullptr passed as pdf\n";
			return 1;
		}
	}
		
	const auto& B_PDFs = Config::getVectorPDFs("Bmass");
        if (int(B_PDFs.size())!=ncontr){
                std::cout<< " NO B_PDFs \n";    
                return 1;
        }
        for (int i=0; i<ncontr; ++i){
                auto pdf = B_PDFs[i].get();
                if (!pdf){
                        std::cout<< "Nullptr passed as pdf\n";
                        return 1;
                }
        }
	std::vector<PDFInterface *> D_PDFs_get;
	std::vector<PDFInterface *> B_PDFs_get;
	for (int i=0; i<ncontr; ++i){
		D_PDFs_get.push_back(D_PDFs[i].get());
		B_PDFs_get.push_back(B_PDFs[i].get());
	}
	//list of fixed variables form config
	for (const auto& fix: fixVect){
		min->FixVariable(min->VariableIndex(fix));
		//std::cout << fix << "  " << min->VariableIndex(fix) << std::endl;
	}
	int randSeed = Config::randSeed;
	TRandom rand;
	if (randSeed >-1)
		rand.SetSeed(randSeed);
	//staring point for stability checks
	for (int ivar=0; ivar<(nvar_md+nvar_mb)*ncontr+ncontr-1; ivar++){
	//for (int ivar=0; ivar<(nvar_md+nvar_mb)*ncontr+ncontr; ivar++){
		double random = 0.0;
		if (randSeed > -1 && !min->IsFixedVariable(ivar))
			random = rand.Uniform(-1.0, 1.0);
		min->SetVariableValue(ivar, starting_point[ivar]*(1.0+0.01*random));
	}

	

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
		replaceIndexVect.push_back(std::make_pair(index_replaced, index_substitute));
		}



	// Start the minimization
	// Define a fit function for Minuit
	auto fchi2 = wrap_chi2(D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB, replaceIndexVect);
	ROOT::Math::Functor f(fchi2, (nvar_md + nvar_mb) * ncontr + ncontr - 1);
	//ROOT::Math::Functor f(fchi2, (nvar_md + nvar_mb) * ncontr + ncontr);
	min->SetFunction(f);
	min->Minimize();

	// Print the fit results
	std::ofstream results(TString::Format("results_%d.txt", sign));
	results << min->Status() << "  " << min->MinValue() << std::endl;
	for (int i = 0; i < (nvar_md + nvar_mb) * ncontr + ncontr - 1; i++)
	//for (int i = 0; i < (nvar_md + nvar_mb) * ncontr + ncontr; i++)
	{
		results << min->X()[i] << "  " << min->Errors()[i] << std::endl;
	}
	results.close();
	return 0;
}
std::function<double(const double*)> wrap_chi2(const std::vector<PDFInterface*>& D_PDFs_get, const std::vector<PDFInterface*>& B_PDFs_get, std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB, const std::vector<std::vector<double>>& dMC_MB, const std::vector<std::pair<int, int>>& replaceIndexVect){

        
	
        const int ncontr = Config::ncontr;

	
	


	double minDM = Config::minDM;
        double maxDM = Config::maxDM;
        double minBMcorr = Config::minBMcorr;
        double maxBMcorr = Config::maxBMcorr;
        const int nvar_md = Config::nvar_md;
        const int nvar_mb = Config::nvar_mb;




	long int emax = vect_2D.size();
	std::cout << emax << " yield of the sample \n";
	auto fchi2 = [D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB, emax, ncontr, minDM, maxDM, minBMcorr, maxBMcorr, nvar_md, nvar_mb, replaceIndexVect](const double *par) -> double 
	{


	// Caclulate chi2
			// Get the pointer in parameters array that correspond to fraction defining parameters????
		const double *pa = &par[ncontr * (nvar_md + nvar_mb)];
		// Calculate fractions
		// Q: Why do we recalculate fractions? What does the minuti minimize - what is stored in `pa` ???
		// The parameters pa[] define the fractions, we have 6 fractions but 5 independent parameters. The parametrisation is arbitrary
		double frac[ncontr];
		double sum_frac = 0.0;
		for (int i=0; i<ncontr-1; i++){
			frac[i] = abs(pa[i]);
			sum_frac +=frac[i];
		}
		frac[ncontr-1] = 1.0-sum_frac;
	
		// Extract the parameters and add some constraints
		double param[ncontr * (nvar_md + nvar_mb)];
		for (int ivar = 0; ivar < (nvar_md+nvar_mb)*ncontr; ivar++)
		{
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var: replaceIndexVect){
			param[irep_var.first] = par[irep_var.second];
		}	
		
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
		if (frac[ncontr -1]>1.0) chi2 += 1.0e6*(frac[ncontr-1]-1.0)*(frac[ncontr-1]-1.0);
		if (frac[ncontr -1]<-0.0) chi2 += 1.0e6*(frac[ncontr-1])*(frac[ncontr-1]);
		return chi2;
	};
	return fchi2;
}

