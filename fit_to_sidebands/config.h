// config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include "json.hpp"

#include "D_M_fit_shape.h"
#include "M_B_2missPT_fit.h"
#include "ChebyshevPDF.h"
#include "pdf_interface.h"

#include <memory>

using namespace cpt_b0_analysis;

//mapping string values of shape vector onto integer values to use switch afterwards
const std::unordered_map<std::string, const int> dictionaryDM ={
	{"DCBplusGaus", 0},
	{"Chebyshev", 1}
};

const std::unordered_map<std::string, const int> dictionaryBMcorr ={
	{"RCplusGaus", 0}
};
const std::unordered_map<std::string, const int> dictionaryChooseFit ={
        {"BM", 0},
	{"DM", 1},
	{"DM+BMfixed", 2},
	{"all", 3}
};


class Config {
public:
	static int load(const std::string& filename);


        static void read_MC(std::vector<std::vector<double>>& xx, std::vector<std::vector<double>>& dxx, std::string MC_directory, int nvar);	
	static std::vector<std::shared_ptr<PDFInterface>> getVectorPDFs(const std::string& domain);

	static bool isMC;
	static bool binned;
	static bool start_from_previous;
    	static int sign;
    	static std::string input_file;
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

    	//Global
    	static double minDM;
    	static double maxDM;
    	static double minBMcorr;
    	static double maxBMcorr;
    	static int nvar_md;
    	static int nvar_mb;
    	static int ncontr;

	static std::vector<std::string> Fits;
	static std::vector<std::string> DMshapes;
	static std::vector<std::string> BMshapes;
	static std::vector <int> int_choose_fits;
	static std::vector <int> intshapesDM;
	static std::vector <int> intshapesBMcorr;

	static std::vector<std::string> fixVect;
	static std::vector<double> fracInit;

	static std::vector<std::string> contrName;
	static std::vector<std::string> varname_md;
	static std::vector<std::string> varname_mb;
	static std::map<std::string, std::string> replace_var;
	static std::map<std::string, std::pair< double, double>> varLimitsMap;

	static std::vector<std::vector<double>> MC_MD;
	static std::vector<std::vector<double>> dMC_MD;
	static std::vector<std::vector<double>> MC_MB;
	static std::vector<std::vector<double>> dMC_MB;

private:
	Config() = default;
};

#endif // CONFIG_H
/* vim:set shiftwidth=8 softtabstop=8 tabstop=8 noexpandtab: */
