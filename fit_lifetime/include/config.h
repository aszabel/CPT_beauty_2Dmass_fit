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
	{"DCBplusGaus", 0}, {"Chebyshev", 1},	{"JSUplusGaus", 2}, {"Exponent", 3},
	{"JSUplusCB", 4},	{"JSUplusDSCB", 5}, {"Johnson", 6},		{"DoubleSidedCrystalBall", 7},
	{"JSUplusExp", 8},	{"DCBplusExp", 9}, {"Gauss", 10}};

const std::unordered_map<std::string, const int> dictionaryBMcorr = {{"RCplusGaus", 0},
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
																	 {"JSUplusDSCB", 14},
																	 {"JSUplusExp", 15},
																	 {"DCBplusExp", 16},
																	 {"Sidebands", 17},
																	 {"Chebyshev", 18}};
const std::unordered_map<std::string, const int> dictionaryChooseFit = {
	{"frac", 0}, {"BM", 1}, {"DM+BMfixed", 2}, {"all", 3}, {"shapes", 4}};

class Config {
public:
	static int load(const std::string& filename);

	static std::vector<std::shared_ptr<PDFInterface>> getVectorPDFs(const std::string& domain);

	static std::string input_file;
	static std::string r_factor_file;

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
	static double r_factor_min;
	static double r_factor_max;
	static double draw_min;
	static double draw_max;

	static double Rez;
	static double Imz;
	static double eta;
	static double alpha;
	static double beta;
	static double t_shift;
	static double imz_start;

	// Global
	static double minDM;
	static double maxDM;
	static double minBMcorr;
	static double maxBMcorr;
	static int Nbins;
	static int nlogBins;
	static std::vector<int> nvar_md;
	static std::vector<int> nvar_mb;
	static std::vector<int> nvar_offset_md;
	static std::vector<int> nvar_offset_mb;
	static int nvar_all_md;
	static int nvar_all_mb;
	static int nvar_time;
	static int ncontr;
	static std::vector<int> ntries;

	static std::vector<std::vector<std::string>> fixVect;

	static std::vector<std::string> contrName;
	static std::vector<std::vector<std::string>> varname_md;
	static std::vector<std::vector<std::string>> varname_mb;
	static std::map<std::string, std::string> replace_var;
	static std::map<std::string, std::pair<double, double>> varLimitsMap;

	static std::vector<std::string> DMshapes;
	static std::vector<std::string> BMshapes;
	static std::vector<int> intshapesDM;
	static std::vector<int> intshapesBMcorr;

private:
	Config() = default;
};

#endif	// CONFIG_H
/* vim:set shiftwidth=8 softtabstop=8 tabstop=8 noexpandtab: */
