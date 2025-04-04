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

#include "TH1D.h"

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
std::function<double(const double*)> wrap_chi2(const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get, const std::vector<std::shared_ptr<PDFInterface>>& B_PDFs_get, std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB, const std::vector<std::vector<double>>& dMC_MB, const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH1D hist1D);



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
        int nentries = Config::nentries;
	

	// Load data set
	TChain ch((Config::chainName).c_str());
	ch.Add(Config::input_file.c_str());
	double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
	double B_MMcorr;
	bool charge;
	int frac_index = 100;
	const int nbins = 50;
	TH1D hist1D("hist1D", "", nbins, Config::minDM, Config::maxDM);

	std::vector<std::pair<double, double>> vect_2D = {};
	if (!Config::isMC){
		ch.SetBranchAddress("B_M", &B_M);
		ch.SetBranchAddress("missPT", &missPT);
		ch.SetBranchAddress("D_M", &D_M);
		ch.SetBranchAddress("mu_PT", &mu_PT);
		ch.SetBranchAddress("mu_P", &mu_P);
		ch.SetBranchAddress("mu_eta", &mu_eta);
		ch.SetBranchAddress("K_PT", &K_PT);
		ch.SetBranchAddress("truecharge", &charge);
	}else{
		ch.SetBranchAddress("BMcorr", &B_MMcorr);
		ch.SetBranchAddress("DM", &D_M);
		ch.SetBranchAddress("frac_index", &frac_index);
	}

	if (nentries < 0) nentries = ch.GetEntries(); 
	if (nentries > ch.GetEntries()){
		std::cerr << "The value of 'nentries' exceeds the number of events in the file." << std::endl;
                return 1;
        }
	int frac_indeces[Config::ncontr];
	for (int i=0; i<Config::ncontr; i++)
		frac_indeces[i] = 0;
	std::cout << nentries << " number of entries read"  <<std::endl;
	for (int i=0; i<nentries; ++i){
		ch.GetEntry(i);
		if (!Config::isMC){
			if (mu_PT < Config::muPTmin || mu_P < Config::muPmin || mu_eta < Config::eta_min || mu_eta > Config::eta_max)
				continue;
			if (int(charge) != Config::sign)
				continue;
			B_MMcorr = B_M + 2.0 * missPT;
		}
		if (B_MMcorr < Config::minBMcorr || B_MMcorr > Config::maxBMcorr)
			continue;
			
		if (D_M < Config::minDM || D_M > Config::maxDM)
			continue;
		vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
		if (D_M<hist1D.GetBinLowEdge(hist1D.FindBin(1830))|| D_M>=hist1D.GetBinLowEdge(hist1D.FindBin(1910)))
			hist1D.Fill(D_M);
		
		if (Config::isMC)
			frac_indeces[frac_index]++;
		
	}


	// Initial fit parameter values taken from 1D fits to MC and Side Bands

		

	// Define Minuit fit variables for M_D
	const int n_all = (Config::nvar_md+Config::nvar_mb) * Config::ncontr + Config::ncontr;
	bool previous_fit = false;
	bool start_scratch = true;
	double starting_point[n_all];
	std::string minName = "Minuit2";
	std::string algoName = "";
	int loop_count = 0;
	for (auto& int_choose_fit: Config::int_choose_fits){
		ROOT::Math::Minimizer *min =
			ROOT::Math::Factory::CreateMinimizer(minName, algoName);

		// Set tolerance , etc...
		min->SetMaxFunctionCalls(Config::functionCalls); // for Minuit/Minuit2
		min->SetTolerance(Config::tolerance[loop_count]);
		loop_count++;
		min->SetPrintLevel(Config::printLevel);
		// min->SetStrategy(2);

		for (int i = 0; i < Config::ncontr; i++)
		{
			for (int ivar = 0; ivar < Config::nvar_md; ivar++)
		{
				min->SetVariable(i * Config::nvar_md + ivar, (TString::Format("%s_%s", Config::contrName[i].c_str(), Config::varname_md[ivar].c_str())).Data(), Config::MC_MD[i][ivar], Config::dMC_MD[i][ivar] + 1.0e-11);
				if (i!=4){
		       			min->FixVariable(i * Config::nvar_md + ivar);
				}else{
					if (ivar>2) 
						min->FixVariable(i * Config::nvar_md + ivar);
				}
				if (start_scratch) 
					starting_point[i * Config::nvar_md + ivar] = Config::MC_MD[i][ivar];

			}
		}

		// Define Minuit fit variables for M_B

/*		for (int i = 0; i < Config::ncontr; i++)
		{
			for (int ivar = 0; ivar < Config::nvar_mb; ivar++)
			{
				min->SetVariable(Config::ncontr * Config::nvar_md + i * Config::nvar_mb + ivar, (TString::Format("%s_%s", Config::contrName[i].c_str(), Config::varname_mb[ivar].c_str())).Data(), Config::MC_MB[i][ivar], Config::dMC_MB[i][ivar] + 1.0e-11);
				if (int_choose_fit == 0 || int_choose_fit == 2){
		       			min->FixVariable(Config::ncontr * Config::nvar_md + i * Config::nvar_mb + ivar);
				}
				if (start_scratch) 
					starting_point[Config::ncontr * Config::nvar_md + i * Config::nvar_mb + ivar] = Config::MC_MB[i][ivar];
			}
		}
*/
		//Set Limits on variables
        	for (auto it =Config::varLimitsMap.begin(); it!=Config::varLimitsMap.end(); ++it)
        	{
                	int index_var = min->VariableIndex(it->first);
                	if (index_var == -1 ){
                        	std::cerr<< "Error in limiting parameters: param " << it->first << " not found.\n";
                        	return 1;
                	}
                	auto pairlims = it->second;
                	min->SetVariableLimits(index_var, pairlims.first, pairlims.second);
        	}

		for (int i = 0; i < Config::ncontr; i++)
		{
			min->SetVariable((Config::nvar_md + Config::nvar_mb) * Config::ncontr + i, (TString::Format("par_frac%d", i)).Data(), Config::fracInit[i], 0.001);
                	//min->SetVariableLimits((Config::nvar_md + Config::nvar_mb) * Config::ncontr + i, -1.0, 1.0);
			//if (int_choose_fit == 1){
			if (i!=4){
				min->SetVariableValue((Config::nvar_md + Config::nvar_mb) * Config::ncontr + i, 0.0);
		       		min->FixVariable((Config::nvar_md + Config::nvar_mb) * Config::ncontr + i);
			}
			if (start_scratch) 
				starting_point[(Config::nvar_md + Config::nvar_mb) * Config::ncontr + i] = Config::fracInit[i];
		}


		// Define the error setimation parameter in minuit for 1 sigma and ncontr -1 free parameters
		double CL_normal = ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1); //1 sigma ~68%
		min->SetErrorDef(TMath::ChisquareQuantile(CL_normal, Config::ncontr-1));// ncontr-1 free fraction parameters, other parameters have gaussian contraints base on MC fits.

		const auto& D_PDFs = Config::getVectorPDFs("Dmass");
		if (int(D_PDFs.size()) != Config::ncontr){
			std::cout<< " NO D_PDFs \n";
			return 1;
		}
		for (int i=0; i<Config::ncontr; ++i){
			auto pdf = D_PDFs[i].get();

			if (!pdf){
				std::cout<< "Nullptr passed as pdf\n";
				return 1;
			}
		}
		
		const auto& B_PDFs = Config::getVectorPDFs("Bmass");
        	/*if (int(B_PDFs.size()) != Config::ncontr){
                	std::cout<< " NO B_PDFs \n";    
                	return 1;
        	}
        	for (int i=0; i<Config::ncontr; ++i){
                	auto pdf = B_PDFs[i].get();
                	if (!pdf){
                        	std::cout<< "Nullptr passed as pdf\n";
                        	return 1;
                	}
        	}
		//list of fixed variables form config*/
		for (const auto& fix: Config::fixVect){
			min->FixVariable(min->VariableIndex(fix));
			//std::cout << fix << "  " << min->VariableIndex(fix) << std::endl;
		}

		/*
		TRandom rand;
		if (Config::randSeed >-1)
			rand.SetSeed(Config::randSeed);
		//staring point for stability checks
		double frac_st[Config::ncontr];
		for (int ivar=0; ivar<(Config::nvar_md+Config::nvar_mb) * Config::ncontr + Config::ncontr; ivar++){
			//for (int ivar=0; ivar<(nvar_md+nvar_mb)*ncontr+ncontr; ivar++){
			double random = 0.0;
			if (Config::randSeed > -1 && !min->IsFixedVariable(ivar) && int_choose_fit==3)
				random = rand.Uniform(-1.0, 1.0);

			min->SetVariableValue(ivar, starting_point[ivar]*(1.0+0.001*random));
		}

		for (int icontr=0; icontr<Config::ncontr; icontr++){
			frac_st[icontr] = starting_point[(Config::nvar_md+Config::nvar_mb) * Config::ncontr + icontr];
			if (Config::randSeed > -1 && !min->IsFixedVariable((Config::nvar_md+Config::nvar_mb) * Config::ncontr + icontr) && int_choose_fit==0){
				double random = rand.Uniform(-1.0, 1.0);
				frac_st[icontr] = random;
			}
		}
			
			//frac_st[0] = 0.5;
			//frac_st[1] = 0.06;
			//frac_st[2] = 0.027;
			//frac_st[3] = 0.06;
			//frac_st[4] = 0.24;
		

		double result0[Config::ncontr*(Config::nvar_md+Config::nvar_mb)+Config::ncontr];
		if (Config::start_from_previous){
			double x, dx;
			std::ifstream res0(Config::previous_result_file.c_str());
			res0>> x >> dx;
			int count = 0;
			while(res0>>x>>dx)
			{
				result0[count] = x;
				count++;
			}
			for (int ic=0; ic<Config::ncontr; ic++){
				double random = 0.0;
				if (Config::randSeed > -1 && !min->IsFixedVariable((Config::nvar_md+Config::nvar_mb) * Config::ncontr+ic) && int_choose_fit==3)
					random = rand.Uniform(-1.0, 1.0);
				frac_st[ic]=result0[(Config::nvar_md+Config::nvar_mb) * Config::ncontr+ic]*(1.0+0.001*random);
			}
		}
		
		double sum_contr=0.0;	
		for (int ic=0; ic<Config::ncontr; ic++)
			sum_contr+= abs(frac_st[ic]);

		for (int ic=0; ic<Config::ncontr; ic++){
			min->SetVariableValue((Config::nvar_md+Config::nvar_mb) * Config::ncontr+ic, frac_st[ic]/sum_contr);
		}
		

	

*/		std::vector<std::pair<int, int>> replaceIndexVect = {};
/*                for (const auto& rep_var: Config::replace_var){
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


*/
		// Start the minimization
		// Define a fit function for Minuit
		auto fchi2 = wrap_chi2(D_PDFs, B_PDFs, vect_2D, Config::MC_MD, Config::dMC_MD, Config::MC_MB, Config::dMC_MB, replaceIndexVect, int_choose_fit, hist1D);
		ROOT::Math::Functor f(fchi2, n_all);
		//ROOT::Math::Functor f(fchi2, (nvar_md + nvar_mb) * ncontr + ncontr);
		min->SetFunction(f);
		min->Minimize();

		// Print the fit results
		std::ofstream results(TString::Format("results1D_DM_%d_%d.txt", Config::sign, int_choose_fit));
		results << min->Status() << "  " << min->MinValue() << std::endl;
		if (min->Status()!=0){
			std::cout << "Bad status of fit "<< min->Status() << std::endl;
			previous_fit = false; 
		}else{
			previous_fit = true;
			start_scratch = false;
		}
		const double *pa = &min->X()[Config::ncontr * (Config::nvar_md + Config::nvar_mb)];
		// Calculate fractions
		// Q: Why do we recalculate fractions? What does the minuti minimize - what is stored in `pa` ???
		// The parameters pa[] define the fractions, we have 6 fractions but 5 independent parameters. The parametrisation is arbitrary
 		double frac_res[Config::ncontr];
                frac_res[0] = abs(pa[0]);
                frac_res[1] = abs(pa[1]);
                frac_res[2] = abs(pa[2]);
                frac_res[3] = abs(pa[3]);
                frac_res[4] = abs(pa[4]);
                frac_res[5] = abs(pa[5]);
			results << min->X()[4 * Config::nvar_md + 0] << "  " << min->Errors()[4 * Config::nvar_md + 0] << std::endl;
			results << min->X()[4 * Config::nvar_md + 1] << "  " << min->Errors()[4 * Config::nvar_md + 1] << std::endl;
			results << min->X()[Config::ncontr * Config::nvar_md + 4] << "  " << min->Errors()[Config::ncontr * Config::nvar_md + 4] << std::endl;

		results.close();
	
		if (min){
			delete min;
		}
		for (int i=0; i<Config::ncontr; i++)
			std::cout << frac_res[i] << "  frac" << i << "  " << double(frac_indeces[i])/double(vect_2D.size()) << std::endl;
	}
	return 0;
}


std::function<double(const double*)> wrap_chi2(const std::vector<std::shared_ptr<PDFInterface>>& D_PDFs_get, const std::vector<std::shared_ptr<PDFInterface>>& B_PDFs_get, std::vector<std::pair<double, double>> vect_2D, const std::vector<std::vector<double>>& MC_MD, const std::vector<std::vector<double>>& dMC_MD, const std::vector<std::vector<double>>& MC_MB, const std::vector<std::vector<double>>& dMC_MB, const std::vector<std::pair<int, int>>& replaceIndexVect, int int_choose_fit, TH1D hist1D){

	long int emax = vect_2D.size();
	std::cout << emax << " yield of the sample \n";
	auto fchi2 = [D_PDFs_get, B_PDFs_get, vect_2D, MC_MD, dMC_MD, MC_MB, dMC_MB, emax, replaceIndexVect, int_choose_fit, hist1D](const double *par) -> double 
	{


	// Caclulate chi2
			// Get the pointer in parameters array that correspond to fraction defining parameters????
		const double *pa = &par[Config::ncontr * (Config::nvar_md + Config::nvar_mb)];
		// Calculate fractions
		// Q: Why do we recalculate fractions? What does the minuti minimize - what is stored in `pa` ???
		// The parameters pa[] define the fractions, we have 6 fractions but 5 independent parameters. The parametrisation is arbitrary
 		double frac[Config::ncontr];

		double chi2 = 0.0;
		//double sum_frac = 0.0;
		for (int i=0; i < Config::ncontr; i++){
			frac[i] = abs(pa[i]);
			//sum_frac+=frac[i];
		}
/*		double tmp = 1.0e5*(sum_frac-1.0)*(sum_frac-1.0);
		if(!Config::binned) 
			tmp *= emax;
		chi2 += tmp;
*/
		//frac[Config::ncontr-1] = abs(1.0-sum_frac);

	        /*frac[0] = 1.0 - abs(pa[0]);
                frac[1] = abs(pa[0]) * (1.0 - abs(pa[1]));
                frac[2] = abs(pa[0]) * abs(pa[1]) * (1.0 - abs(pa[2]));
                frac[3] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * (1.0 - abs(pa[3]));
                frac[4] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * abs(pa[3]) * (1.0 - abs(pa[4]));
                frac[5] = abs(pa[0]) * abs(pa[1]) * abs(pa[2]) * abs(pa[3]) * abs(pa[4]);
	*/
		// Extract the parameters and add some constraints
		double param[Config::ncontr * (Config::nvar_md + Config::nvar_mb)];
		for (int ivar = 0; ivar < (Config::nvar_md+Config::nvar_mb)*Config::ncontr; ivar++)
		{
			param[ivar] = par[ivar];
		}
		for (const auto& irep_var: replaceIndexVect){
			param[irep_var.first] = par[irep_var.second];
		}	
		
		// Calculate normalisation integrals
		for (int i = 0; i < Config::ncontr; i++)
		{
			D_PDFs_get[i]->CalcIntegral(&param[i * Config::nvar_md], Config::minDM, Config::maxDM);
			//B_PDFs_get[i]->CalcIntegral(&param[Config::ncontr * Config::nvar_md + i * Config::nvar_mb], Config::minBMcorr, Config::maxBMcorr);
		}
/*
		if(int_choose_fit != 1){ //don't calculate if only Bmass fit
			for (int i = 0; i < Config::ncontr; i++)
			{
				// TODO set dxx for D comb background to 0 and remove this
				for (int ivar = 0; ivar < Config::nvar_md; ivar++)
				{
					// Skip the constrain for chebyshev background: {"Chebyshev", 1}
					if (Config::intshapesDM[i] == 1 || ivar == 1) continue;
					if (dMC_MD[i][ivar] != 0){
					 	double tmp = (MC_MD[i][ivar] - param[i * Config::nvar_md + ivar]) * (MC_MD[i][ivar] - param[i * Config::nvar_md + ivar]) / (2.0*(dMC_MD[i][ivar] * dMC_MD[i][ivar])); // use the results of MC fits
						chi2 += 2.5e3*tmp;
					}
				}
			}
		}
		if (int_choose_fit != 1){ //don't calculate if only Dmass fit
			for (int i = 0; i < Config::ncontr; i++)
		{
				for (int ivar = 0; ivar < Config::nvar_mb; ivar++)
				{
					if (dMC_MB[i][ivar] != 0){
						double tmp = (MC_MB[i][ivar] - param[Config::ncontr * Config::nvar_md + i * Config::nvar_mb + ivar]) * (MC_MB[i][ivar] - param[Config::ncontr * Config::nvar_md + i * Config::nvar_mb + ivar]) /(2.0* (dMC_MB[i][ivar] * dMC_MB[i][ivar])); // use the results of MC fits																							    		    	

						chi2 += 2.5e3*tmp;
					}
				}
			}
		}
*/
		double *vect_chi2 = new double[vect_2D.size()]; // Per event results - required to efficiently calculate a Kahan compensated sum
		for (int jj = 0; jj<int(vect_2D.size()); jj++)
				vect_chi2[jj] = 0.0;
		double test_chi2 = 0.0; 
		if (Config::binned)
		{
			int nbins_md = hist1D.GetNbinsX();
			//int nbins_mb = 0;
			//double bin_width_mb = (Config::maxBMcorr-Config::minBMcorr)/double(nbins_mb);
			double bin_width_md = (Config::maxDM-Config::minDM)/double(nbins_md);
			//double histev = hist1D.Integral();
#pragma omp parallel for
			for (int bin_md=1; bin_md<=nbins_md; bin_md++)
			{
				//for (int bin_mb=1; bin_mb<=nbins_mb; bin_mb++)
				//{
					double sum_contr = 0.0;
					for (int i = 0; i < Config::ncontr; i++)
                        		{
						double mdass = hist1D.GetXaxis()->GetBinCenter(bin_md);
						//double mcorr = hist1D.GetYaxis()->GetBinCenter(bin_mb);
						double md_val = D_PDFs_get[i]->EvalPDF(&mdass, &param[i * Config::nvar_md]);
						//double mb_val = B_PDFs_get[i]->EvalPDF(&mcorr, &param[Config::ncontr * Config::nvar_md + i * Config::nvar_mb]);
						sum_contr += emax*bin_width_md*frac[i]*md_val; 
						//std::cout << sum_contr << "  sumcontr  " << frac[i] << "  " << md_val << "  " << mb_val << std::endl;
					}
					double cont = (double)hist1D.GetBinContent(bin_md);
					double err  = (double)hist1D.GetBinError(bin_md);
					if (err !=0.0)
					{ 
						vect_chi2[bin_md-1]+=0.5*(sum_contr-cont)*(sum_contr-cont)/err/err;
						test_chi2+=0.5*(sum_contr-cont)*(sum_contr-cont)/err/err;
						//std::cout << bin_md << "   "  << bin_mb << "  " << cont << "  " << sum_contr << "  " << sum_contr/cont << std::endl; 
					}
					//std::cout << test_chi2 << " chi2 \n";
				//}
			}
		}else{
					
		// Main loop that calculates the chi2
// Main loop that calculates the chi2 - run in parallel using OpenMP
#pragma omp parallel for
		for (long int e = 0; e < emax; e++)
		{
			double mdass = std::get<0>(vect_2D[e]);
			double mcorr = std::get<1>(vect_2D[e]);
			double like_event = 0.0;
			for (int i = 0; i < Config::ncontr; i++)
			{
				double md_like, mb_like;
				md_like = D_PDFs_get[i]->EvalPDF(&mdass, &param[i * Config::nvar_md]);
				mb_like = B_PDFs_get[i]->EvalPDF(&mcorr, &param[Config::ncontr * Config::nvar_md + i * Config::nvar_mb]);
				if (mb_like <0.0 || mb_like>1.0 || md_like <0.0 || md_like>1.0 || frac[i] <0.0 || frac[i] >1.0){
					chi2+=1e15;
					continue;
				}
				switch(int_choose_fit){
					case 0:
						like_event += mb_like * frac[i];
						break;
					case 1:
						like_event += md_like * frac[i];
						break;
				        case 2:	
						like_event += md_like * mb_like * frac[i];
						break;
				        case 3:	
						like_event += md_like * mb_like * frac[i];
						break;
					default:
						std::cerr << "Wrong fit choice!" << std::endl;
					       	break;
				}
				
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
		if (avx)
		{
			size_t vecsize;
			if (Config::binned)
				vecsize = size_t(hist1D.GetNbinsX()); 
			else
				vecsize = (size_t)vect_2D.size();
			chi2_threads = fastAccurate<Method::Kahan, 4>(vect_chi2, vecsize);
		}
		else
		{
			double sum = 0, c = 0;
			for (long unsigned int i = 0; i < vect_2D.size(); i++)
			{
				ksum(sum, c, vect_chi2[i]);
			}
			chi2_threads = sum + c;
		}
		delete [] vect_chi2;
		chi2 += chi2_threads;
		// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		// std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
		if (!Config::binned) chi2 += log(double(emax));		//Extended ML 

		 
		
		return chi2;
	};
	return fchi2;
}

/* vim:set shiftwidth=8 softtabstop=8 tabstop=8 noexpandtab: */
