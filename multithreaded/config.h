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

class Config {
public:
	static int load(const std::string& filename);


        static void read_MC(std::vector<std::vector<double>>& xx, std::vector<std::vector<double>>& dxx, std::string MC_directory, int nvar);	
	static std::vector<std::unique_ptr<PDFInterface>> getVectorPDFs(const std::string& domain);

    	static int sign;
    	static std::string input_file;

    	static int nentries;

    	static double tolerance;

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
    	static int nbins;

	static std::vector<std::string> DMshapes;
	static std::vector<std::string> BMshapes;

	static std::vector<std::string> fixVect;
	static std::vector<double> fracInit;

	static std::vector<std::string> contrName;
	static std::vector<std::string> varname_md;
	static std::vector<std::string> varname_mb;

	static std::vector<std::vector<double>> MC_MD;
	static std::vector<std::vector<double>> dMC_MD;
	static std::vector<std::vector<double>> MC_MB;
	static std::vector<std::vector<double>> dMC_MB;

private:
    	Config() = default;
	static std::vector <int> intshapesDM;
        static std::vector <int> intshapesBMcorr;

};

#endif // CONFIG_H
