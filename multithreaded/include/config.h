// config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <fstream>
#include <memory>
#include <string>

#include "BasicShapes.h"
#include "ChebyshevPDF.h"
#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "json.hpp"
#include "pdf_interface.h"

using namespace cpt_b0_analysis;

// mapping string values of shape vector onto integer values to use switch
// afterwards
const std::unordered_map<std::string, const int> dictionaryDM = {
	{"DCBplusGaus", 0}, {"Chebyshev", 1}, {"JSUplusGaus", 2},
	{"Exponent", 3},	{"JSUplusCB", 4}, {"JSUplusDSCB", 5},
	{"Johnson", 6},
	{"DoubleSidedCrystalBall", 7}
};

const std::unordered_map<std::string, const int> dictionaryBMcorr = {
	{"RCplusGaus", 0},
	{"SNplusCB", 1},
	{"Gauss", 2},
	{"DoubleGauss", 3},
	{"CrystalBall", 4},
	{"SkewNormal", 5},
	{"RaisedCosine", 6},
	{"Johnson", 7},
	{"DoubleSidedCrystalBall", 8},
	{"JSUplusCB", 9},
	{"DCBplusGaus", 10},
	{"JSUplusGaus", 11},
	{"JSUplusCBplusDGaus", 12},
	{"SNplusCBplusDGaus", 13},
};
const std::unordered_map<std::string, const int> dictionaryChooseCategory = {
	{"2D", 0},
	{"1D_DM", 1},
	{"1D_BM", 2},
};
const std::unordered_map<std::string, const int> dictionaryChooseMassVariable = {
	{"D_M", 0}, {"B_M", 1}, {"B_Mcorr", 2}, {"B_MMcorr", 3}, {"missPT", 4}};
const std::unordered_map<std::string, const int> dictionaryChooseFit = {
	{"frac", 0}, {"BM", 1}, {"DM+BMfixed", 2}, {"all", 3}, {"shapes", 4}};
const std::unordered_map<std::string, const int> dictionaryChooseFit1D = {
	{"signal", 0},	 {"BuDmunu", 1},   {"BsDsMunu", 2},
	{"B02DpDsm", 3}, {"sidebands", 4}, {"Bu2D0Dsm", 5},
};

class Config {
public:
	static int load(const std::string& filename);

	static void read_MC(std::vector<std::vector<double>>& xx, std::vector<std::vector<double>>& dxx,
						std::string MC_directory, std::vector<int> nvar);
	static std::vector<std::shared_ptr<PDFInterface>> getVectorPDFs(const std::string& domain);

	static bool isMC;
	static bool is2D;
	static bool binned;
	static bool start_from_previous;
	static bool calc_sWeights;
	static int sign;
	static std::string category;
	static int int_category;
	static std::string mass_variable;
	static int int_mass_variable;
	static std::vector<std::string> input_files;
	static std::string previous_result_file;

	static int nentries;

	static std::vector<double> tolerance;
	static int functionCalls;
	static int printLevel;
	static int randSeed;

	static std::string chainName;
	static double muPTmin;
	static double muPmin;
	static double eta_min;
	static double eta_max;
	static double tMin;
	static double tMax;

	// Global
	static double minDM;
	static double maxDM;
	static double minBMcorr;
	static double maxBMcorr;
	static std::vector<int> nvar_md;
	static std::vector<int> nvar_mb;
	static std::vector<int> nvar_offset_md;
	static std::vector<int> nvar_offset_mb;
	static int nvar_all_md;
	static int nvar_all_mb;
	static int ncontr;
	static int n_sideband;
	static int ntries;

	static std::vector<std::string> Fits;
	static std::vector<std::string> DMshapes;
	static std::vector<std::string> BMshapes;
	static std::vector<int> int_choose_fits;
	static std::vector<int> intshapesDM;
	static std::vector<int> intshapesBMcorr;
	static std::vector<int> test;

	static std::vector<std::string> fixVect;
	static std::vector<double> fracInit;
	static std::vector<std::vector<double>> init_values;

	static std::vector<std::string> contrName;
	static std::vector<std::vector<std::string>> varname_md;
	static std::vector<std::vector<std::string>> varname_mb;
	static std::map<std::string, std::string> replace_var;
	static std::map<std::string, std::pair<double, double>> varLimitsMap;

	static std::vector<std::vector<double>> MC_MD;
	static std::vector<std::vector<double>> dMC_MD;
	static std::vector<std::vector<double>> MC_MB;
	static std::vector<std::vector<double>> dMC_MB;

	static std::string MC_directory_MD;
	static std::string MC_directory_MB;

private:
	Config() = default;
};

#endif	// CONFIG_H
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
