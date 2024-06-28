#include <string>
#include <fstream>
#include <iomanip>
#include "omp.h"

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

// TODO is this realy required??
typedef std::basic_string<char> string;
typedef std::basic_ifstream<char> ifstream;
typedef std::basic_ofstream<char> ofstream;

using namespace cpt_b0_analysis;

#include <mutex>
std::mutex my_mutex;

// TODO move those globals to a config file
double minx = 1800.;
double maxx = 1940.;
double miny = 2700.;
double maxy = 8300.;
const int nvar_md = 7;
const int nvar_mb = 7;
const int ncontr = 6;
const int nbins = 40;
bool uncorr = false;
//bool avx = false;
bool avx = true;

/**
 * @brief Kahan compensated summation based on http://blog.zachbjornson.com/2019/08/11/fast-float-summation.html
 * 
 * @param sum Load and store value of the kahan sum
 * @param c Load and store value of the summation error
 * @param y Value to add
 */
inline static void kadd(double& sum, double& c, double y) {
  y -= c;
  auto t = sum + y;
  c = (t - sum) - y;
  sum = t;
}

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
int main(int argc, char *argv[])
{
	// The first, fundamental operation to be performed in order to make ROOT
	// thread-aware.
	ROOT::EnableThreadSafety();

	// Increase printout precision
	//std::cout<<std::setprecision(40);

	if (argc != 2)
	{
		std::cout << " Please state: 'muplus' or 'muminus'." << std::endl;
		return 1;
	}
	if (strcmp(argv[1], "muplus") != 0 && strcmp(argv[1], "muminus") != 0)
	{
		std::cout << " Please state: 'muplus' or 'muminus'." << std::endl;
		return 1;
	}
	int sign;
	if (strcmp(argv[1], "muplus") == 0)
		sign = 1;
	if (strcmp(argv[1], "muminus") == 0)
		sign = 0;

	string minName = "Minuit2";
	string algoName = "";
	ROOT::Math::Minimizer *min =
		ROOT::Math::Factory::CreateMinimizer(minName, algoName);

	// Set tolerance , etc...
	// TODO load those values from a config file
	min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
	min->SetMaxIterations(10000);		// for GSL
	// min->SetTolerance(2.50e5);
	min->SetTolerance(50.);
	min->SetPrintLevel(2);
	// min->SetStrategy(2);
	// min->SetPrecision(0.00001);

	// Load data set
	TChain ch("BlindedTree");
	// TODO load the data filename from config file
	ch.Add("/mnt/home/share/lhcb/CPT_beauty/data2016/selected/selected_data2016MagDown.root");
	// ch.Add("/home/szabelskia/LHCb/data2016/tree_missPT_D_M_MagDown25102023_nomassDmuCut/selected_data2016MagDown.root");
	double D_M, mu_PT, mu_P, mu_eta, K_PT, B_M, missPT;
	bool charge;

	std::vector<std::pair<double, double>> vect_2D;
	TH2D *hist = new TH2D("hist", "", nbins, minx, maxx, nbins, miny, maxy);
	ch.SetBranchAddress("B_M", &B_M);
	ch.SetBranchAddress("missPT", &missPT);
	ch.SetBranchAddress("D_M", &D_M);
	ch.SetBranchAddress("mu_PT", &mu_PT);
	ch.SetBranchAddress("mu_P", &mu_P);
	ch.SetBranchAddress("mu_eta", &mu_eta);
	ch.SetBranchAddress("K_PT", &K_PT);
	ch.SetBranchAddress("truecharge", &charge);

	for (int i=0; i<ch.GetEntries(); ++i)
	// TODO load number of entries from the config file
	//for (int i = 0; i < 1e5; ++i)
	{
		ch.GetEntry(i);
		// TODO load cut values from a config file
		if (mu_PT < 500 || mu_P < 5000 || mu_eta < 2 || mu_eta > 4.5)
			continue;
		double B_MMcorr = B_M + 2.0 * missPT;
		if (B_MMcorr < miny || B_MMcorr > maxy)
			continue;
		if (D_M < minx || D_M > maxx)
			continue;
		if (int(charge) == sign)
		{
			vect_2D.push_back(std::make_pair(D_M, B_MMcorr));
			hist->Fill(D_M, B_MMcorr);
		}
	}
	// double  nevents = double(vect_2D.size());

	// TODO define shapes dynamically preferably based on configfile??
	md_fit md_shape(minx, maxx, ncontr);
	mb_2misspt_fit mb_shape(miny, maxy, ncontr);

	// Initial fit parameter values taken from 1D fits to MC and Side Bands
	double xx[ncontr][nvar_md];
	double dxx[ncontr][nvar_md];
	// TODO define this in a header with some global params that describe the fit??
	TString name[ncontr] = {"signal", "BuDmunu", "BsDsMunu", "B02DpDsm", "sidebands", "Bu2D0Dsm"};
	// Read results of 2D fits that are stored as simple text files
	// Each line corresponds to a single parameter
	// First column defines parameter value
	// Second column defines parameter uncertainty
	// D mass fit
	for (int ifile = 0; ifile < ncontr; ifile++)
	{
		std::ifstream infile(TString::Format("D_M_results/res_%s_%d.txt", name[ifile].Data(), sign));
		int k = 0;
		while (infile >> xx[ifile][k] >> dxx[ifile][k])
		{
			k++;
		}
	}
	
	// B mass fit
	double xx_mcorr[ncontr][nvar_mb];
	double dxx_mcorr[ncontr][nvar_mb];
	for (int ifile = 0; ifile < ncontr; ifile++)
	{
		ifstream infile_mb(TString::Format("B_M_results/res_%s_%d.txt", name[ifile].Data(), sign));
		int k = 0;
		while (infile_mb >> xx_mcorr[ifile][k] >> dxx_mcorr[ifile][k])
		{
			k++;
		}
	}

	// Caclulate chi2
	// TODO consider changing this to a full function
	auto fchi2 = [&md_shape, &mb_shape, vect_2D, xx, dxx, xx_mcorr, dxx_mcorr](const double *par) -> double
	{
		// Get the pointer in parameters array that corresponds to location just after shape parameters - number of events ????
		const double *pa = &par[ncontr * (nvar_md + nvar_mb)];
		// Calculate fractions
		// Q: Why do we recalculate fractions? What does the minuti minimize - what is stored in `pa` ???
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
				if (ivar == 1 && i != 2 && i != 4)
					param[i * nvar_md + ivar] = par[1];
			}
		}
		// Q: Second set of params ???? Number of events ???
		for (int i = 0; i < ncontr; i++)
		{
			for (int ivar = 0; ivar < nvar_mb; ivar++)
			{
				int ipar = ncontr * nvar_md + i * nvar_mb + ivar;
				param[ipar] = par[ipar];
			}
		}
		/*
		for (int i=0; i<ncontr * (nvar_md + nvar_mb); i++) {
			std::cout<<param[i]<<", ";
		}
		std::cout<<pa[0]<<", "<<pa[1]<<", "<<pa[2]<<", "<<pa[3]<<std::endl;
		*/

		// Calculate normalisation integrals
		for (int i = 0; i < ncontr; i++)
		{
			// TODO move to proper class
			double sigma = param[i * nvar_md];
			double mean = param[i * nvar_md + 1];
			md_shape.int_gaus.push_back(ROOT::Math::normal_cdf(maxx, sigma, mean) - ROOT::Math::normal_cdf(minx, sigma, mean));
			sigma = param[i * nvar_md + 2];
			double n = param[i * nvar_md + 5];
			double alpha = param[i * nvar_md + 4];
			double alpha_h = param[i * nvar_md + 6];

			md_shape.int_DCB[i] = TMath::Abs(-ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean) + ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2. * mean - maxx, alpha_h, n, sigma, mean) + ROOT::Math::crystalball_integral(mean, alpha_h, n, sigma, mean));

			if (i == 4)
			{
				double a1 = par[i * nvar_md + 0];
				double a2 = par[i * nvar_md + 1];

				md_shape.int_Cheb = abs((1.0 - a2) * maxx + 0.5 * a1 * maxx * maxx + 2. / 3. * a2 * maxx * maxx * maxx - (1.0 - a2) * minx - 0.5 * a1 * minx * minx - 2. / 3. * a2 * minx * minx * minx);
			}

			double s = abs(param[ncontr * nvar_md + i * nvar_mb + 1]);
			mean = param[ncontr * nvar_md + i * nvar_mb];
			double minnB = mean - s;
			if (minx > mean - s)
				minnB = miny;
			double maxxB = mean + s;
			if (maxx < mean + s)
				maxxB = maxy;

			mb_shape.int_cos[i] = 1.0 / (2.0) * (1.0 + (maxxB - mean) / s + TMath::Sin((maxxB - mean) / s * TMath::Pi()) / TMath::Pi()) - 1.0 / (2.0) * (1.0 + (minnB - mean) / s + TMath::Sin((minnB - mean) / s * TMath::Pi()) / TMath::Pi());

			double mean1 = param[ncontr * nvar_md + i * nvar_mb + 3];
			double sigma1 = abs(param[ncontr * nvar_md + i * nvar_mb + 4]);

			mb_shape.int_gaus1[i] = ROOT::Math::normal_cdf(maxy, sigma1, mean1) - ROOT::Math::normal_cdf(miny, sigma1, mean1);

			double mean2 = param[ncontr * nvar_md + i * nvar_mb + 5];
			double sigma2 = abs(param[ncontr * nvar_md + i * nvar_mb + 6]);

			mb_shape.int_gaus2[i] = ROOT::Math::normal_cdf(maxy, sigma2, mean2) - ROOT::Math::normal_cdf(miny, sigma2, mean2);
		}

		// Main loop that calculates the chi2
		double chi2 = 0.0;  // Final result		
		double* vect_chi2 = new double[vect_2D.size()];  // Per event results - required to efficiently calculate a Kahan compensated sum

		// Main look - run in parallel using OpenMP
		#pragma omp parallel for
		for (long unsigned int i=0; i<vect_2D.size(); i++)
		{
			/*
			// Use mutex to have sensible printouts from the parallel section
			{
				std::lock_guard<std::mutex> guard0(my_mutex);
				std::cout<<"Chi2 before: "<<chi2_tmp<<std::flush<<std::endl;
			}
			*/
			double mdass = std::get<0>(vect_2D[i]);
			double mcorr = std::get<1>(vect_2D[i]);
			double likelihood = 0.0;
			for (int k = 0; k < ncontr; k++)
			{
				double md_like, mb_like;
				md_like = md_shape.func_full(&mdass, &param[k * nvar_md], k);
				mb_like = mb_shape.func_full(&mcorr, &param[ncontr * nvar_md + k * nvar_mb], k);
				likelihood += md_like * mb_like * frac[k];
			}
			if (likelihood > 0.0)
				vect_chi2[i] = - 2.0 * log(likelihood);
			/*
			{
				std::lock_guard<std::mutex> guard(my_mutex);
				std::cout<<"Chi2: "<<chi2_tmp<<std::flush<<std::endl;
			}
			*/
		}

		//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		if (avx) {
			chi2 = fastAccurate<Method::Kahan, 4>(vect_chi2, (size_t)vect_2D.size());
		} else {
			double sum = 0, c = 0;
  			for (long unsigned int i=0; i<vect_2D.size(); i++) {
    			kadd(sum, c, vect_chi2[i]);
  			}
  			chi2 = sum + c;
		}
		delete vect_chi2;
		//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;

		// Add additional constraints based on the templates from 1D fits
		for (int i = 0; i < ncontr; i++)
		{
			if (i == 4)
				continue;
			for (int ivar = 0; ivar < nvar_md; ivar++)
			{
				if (dxx[i][ivar] != 0)
					chi2 += (xx[i][ivar] - param[i * nvar_md + ivar]) * (xx[i][ivar] - param[i * nvar_md + ivar]) / (dxx[i][ivar] * dxx[i][ivar]); // use the results of MC fits
			}
		}
		for (int i = 0; i < ncontr; i++)
		{
			for (int ivar = 0; ivar < nvar_mb; ivar++)
			{
				if (dxx_mcorr[i][ivar] != 0)
					chi2 += (xx_mcorr[i][ivar] - param[ncontr * nvar_md + i * nvar_mb + ivar]) * (xx_mcorr[i][ivar] - param[ncontr * nvar_md + i * nvar_mb + ivar]) / (dxx_mcorr[i][ivar] * dxx_mcorr[i][ivar]); // use the results of MC fits
			}
		}

		//std::cerr<<"Chi2: "<<chi2<<std::endl;
		return chi2;
	};

	// Define Minuit fit variables for M_D
	TString varname_md[nvar_md] = {"sigma", "mean", "sigmaCB", "f12", "alpha", "n", "alpha_h"};

	for (int i = 0; i < ncontr; i++)
	{
		for (int ivar = 0; ivar < nvar_md; ivar++)
		{
			min->SetVariable(i * nvar_md + ivar, (TString::Format("%s_%s", name[i].Data(), varname_md[ivar].Data())).Data(), xx[i][ivar], dxx[i][ivar] + 1.0e-11);
			// min->FixVariable(i*nvar_md+ivar);
			min->SetVariableLowerLimit(i * nvar_md + ivar, 1e-13);
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

	// Define Minuit fit variables for M_D
	TString varname_mb[nvar_mb] = {"mean_rc", "sigma_rc", "f12_gaus", "mean1_gaus", "sigma1_gaus", "mean2_gaus", "sigma2_gaus"};

	for (int i = 0; i < ncontr; i++)
	{
		// TODO load limits and fixed value from the config file
		for (int ivar = 0; ivar < nvar_mb; ivar++)
		{
			min->SetVariable(ncontr * nvar_md + i * nvar_mb + ivar, (TString::Format("%s_%s", name[i].Data(), varname_mb[ivar].Data())).Data(), xx_mcorr[i][ivar], dxx_mcorr[i][ivar] + 1.0e-11);
			// min->FixVariable(ncontr*nvar_md+i*nvar_mb+ivar);
			min->SetVariableLowerLimit(ncontr * nvar_md + i * nvar_mb + ivar, 1e-13);
			if (ivar == 2)
				min->SetVariableLimits(ncontr * nvar_md + i * nvar_mb + ivar, -1.0, 1.0);
		}
	}

	// Set initial fraction values
	// TODO load this from a config file
	double frac_init[ncontr - 1] = {0.468636, 0.736568, 0.973902, 0.949009, -0.00476566};
	double frac_init_plus[ncontr - 1] = {0.578309, 0.816153, 0.519927, 0.627151, 0.0151929};

	for (int i = 0; i < ncontr - 1; i++)
	{

		if (sign == 0)
			min->SetVariable((nvar_md + nvar_mb) * ncontr + i, (TString::Format("par_frac%d", i)).Data(), frac_init[i], 0.01);
		else
			min->SetVariable((nvar_md + nvar_mb) * ncontr + i, (TString::Format("par_frac%d", i)).Data(), frac_init_plus[i], 0.01);

		min->SetVariableLimits((nvar_md + nvar_mb) * ncontr + i, -1.0, 1.0);
	}

	// Define a fit function for Minuit
	ROOT::Math::Functor f(fchi2, (nvar_md + nvar_mb) * ncontr + ncontr - 1);
	min->SetFunction(f);

	// Q: Define ??????
	double CL_normal = ROOT::Math::normal_cdf(1) - ROOT::Math::normal_cdf(-1);
	min->SetErrorDef(TMath::ChisquareQuantile(CL_normal, min->NFree()));

	// Q: Load values of fixed parameters ?????
	double x_fix[ncontr * (nvar_md + nvar_mb) + ncontr - 1];
	double dx_fix[ncontr * (nvar_md + nvar_mb) + ncontr - 1];
	ifstream input_fix(TString::Format("results_fix_%d.txt", sign));
	int k = 0;
	while (input_fix >> x_fix[k] >> dx_fix[k])
		k++;
	input_fix.close();

	// Q: WHAT IS THIS ?? What is the purpose of the `uncorr` flag?
	uncorr = true;
	if (k == ncontr * (nvar_md + nvar_mb) + ncontr - 1)
		std::cout << "TUUUUUUU \n\n";

	/*
	for(int ivar=0; ivar<ncontr*(nvar_md+nvar_mb)+ncontr-1; ivar++){
		std::cout << x_fix[ivar] << "  " << ivar << std::endl;
		min->SetVariableValue(ivar, x_fix[ivar]);
	}
	*/

	// Start the minimization
	min->Minimize();

	// Verify the parameter correlations
	// Q: WHY IS `over` not used ????
	bool over[ncontr * (nvar_md + nvar_mb)];
	for (int i = 0; i < ncontr * (nvar_md + nvar_mb); i++)
	{
		over[i] = false;
		for (int j = 0; j < ncontr - 1; j++)
		{
			double cor = min->Correlation(i, ncontr * (nvar_md + nvar_mb) + j);
			if (abs(cor) < 0.02)
			{
				std::cout << min->VariableName(i) << "  " << cor << std::endl;
			}
			else
				over[i] = true;
		}
	}

	// Q: Fix parameters for the sFit ?????
	// TODO think of a way to define it separately via config file or a header
	// Q: What do we fix here ?? 
	std::vector<int> vect_fix = {5, 8, 12, 15, 19, 22, 26, 29, 30, 31, 32, 33, 34, 36, 40, 57, 62};
	for (int i = 0; i < ncontr * (nvar_md + nvar_mb) + ncontr - 1; i++)
		min->SetVariableValue(i, x_fix[i]);
	uncorr = true;
	for (auto i : vect_fix)
		if (uncorr)
			min->FixVariable(i);

	// min->Minimize();

	// Print the fit results
	ofstream results(TString::Format("results_%d.txt", sign));
	for (int i = 0; i < (nvar_md + nvar_mb) * ncontr + ncontr - 1; i++)
	{
		results << min->X()[i] << "  " << min->Errors()[i] << std::endl;
	}
	results.close();

	return 0;
}
