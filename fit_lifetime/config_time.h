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

        const std::unordered_map<std::string, const int> dictionaryDM ={
                {"DCBplusGaus", 0},
                {"Chebyshev", 1}
        };

        const std::unordered_map<std::string, const int> dictionaryBMcorr ={
                {"RCplusGaus", 0}
        };
        const std::unordered_map<std::string, const int> dictionaryChooseFit ={
                {"frac", 0},
                {"BM", 1},
                {"DM+BMfixed", 2},
                {"all", 3},
                {"shapes", 4}
        };



class Config {
public:
	static int load(const std::string& filename);


	static std::vector<std::shared_ptr<PDFInterface>> getVectorPDFs(const std::string& domain);

    	static std::string input_file;
    	static std::string r_factor_file;


    	static int nentries;

    	static double tolerance;
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

    	//Global
    	static double minDM;
    	static double maxDM;
    	static double minBMcorr;
    	static double maxBMcorr;
    	static int Nbins;
    	static int nlogBins;
    	static int nvar_md;
    	static int nvar_mb;
    	static int nvar_time;
    	static int ncontr;
    	static int ntries;


	static std::vector<std::vector<std::string>> fixVect;

	static std::vector<std::pair<std::string, double>> varname;
	static std::map<std::string, std::string> replace_var;
	static std::map<std::string, std::pair< double, double>> varLimitsMap;

        static std::vector<std::string> DMshapes;
        static std::vector<std::string> BMshapes;
        static std::vector <int> intshapesDM;
        static std::vector <int> intshapesBMcorr;

private:
	Config() = default;
};

#endif // CONFIG_H
/* vim:set shiftwidth=8 softtabstop=8 tabstop=8 noexpandtab: */
