// config.cpp
#include "config.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using json = nlohmann::json;
using namespace cpt_b0_analysis;


int Config::sign = -1;
std::string Config::input_file = "";

int Config::nentries = -1;

double Config::tolerance = 50.;

double Config::muPTmin = 500.;
double Config::muPmin = 5000.;
double Config::eta_min  = 2.;
double Config::eta_max = 4.5;

//Global
double Config::minDM = 1800.;
double Config::maxDM = 1940.;
double Config::minBMcorr = 2700.;
double Config::maxBMcorr = 8300.;
int Config::nvar_md = 7;
int Config::nvar_mb = 7;
int Config::ncontr = 6;
int Config::nbins = 40;

std::vector<int> Config::intshapesDM = {};
std::vector<int> Config::intshapesBMcorr = {};

std::vector<std::string> Config::DMshapes = {};
std::vector<std::string> Config::BMshapes = {};

std::vector<std::string> Config::fixVect = {};

std::vector<std::unique_ptr<PDFInterface>> Config::getVectorPDFs(const std::string& domain) {
	std::vector<std::unique_ptr<PDFInterface>> vPDFs;
	std::cout << domain << std::endl;
	if (domain == std::string("Dmass")){
		std::cout << intshapesDM.size() << " size\n";
		for (auto intshape: intshapesDM){
				std::cout << intshape << "TUTAJ" << std::endl;
                	switch (intshape){
                        	case 0:
                                	vPDFs.push_back(std::make_unique<DoubleSidedCrystalballPlusGaussPDF>());
                                	break;
                        	case 1:
                                	vPDFs.push_back(std::make_unique<ChebyshevPDF>());
                                	break;
                        	default:
                                	std::cerr << "Error while DM_pdf dynamic declaration. Check if shapes from config file refer to the shapes defined in the code." << std::endl;
                                	return std::vector<std::unique_ptr<PDFInterface>>{};
                	}
		}

	}else if (domain == std::string("Bmass")){
	        for (auto intshape: intshapesBMcorr){
			std::cout << intshape << " tam " <<std::endl;
        	        switch (intshape){
                        	case 0:
                                	vPDFs.push_back(std::make_unique<RaisedCosinePlusGaussPDF>());
                                	break;
                        	default:
                                	std::cerr << "Error while BM_pdf dynamic declaration. Check if shapes from config file refer to the shapes defined in the code." << std::endl;
                                return std::vector<std::unique_ptr<PDFInterface>>{};
                	}
        	}

	}else{
		std::cerr << "Error wrong domaine in getVectorPDFs(domain) chose 'Dmass' or 'Bmass'" << std::endl;
		return std::vector<std::unique_ptr<PDFInterface>>{};
	}	
    	return vPDFs;
}

int Config::load(const std::string& filename) {
	std::ifstream config_file(filename);
        json config = json::parse(config_file);
        config_file.close();

        if (config.contains("sign")) {
                if (config["sign"].template get<std::string>() == std::string("muplus"))
                        sign = 1;
                else if (config["sign"].template get<std::string>() == std::string("muminus"))
                        sign = 0;
                else {
                        std::cout << "Invalid config file: allowed values for the 'sign' key are 'muplus' or 'muminus'." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'sign' key." << std::endl;
                return 1;
        }


        if (config.contains("tolerance")){
                tolerance = config["tolerance"];
                if (tolerance <0.){
                        std::cout << "Invalid config file: 'tolerance' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'tolerance' key." << std::endl;
                return 1;
        }


        if (config.contains("muPmin")){
                muPmin = config["muPmin"];
                if (muPmin <0){
                        std::cout << "Invalid config file: 'muPmin' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'muPmin' key." << std::endl;
                return 1;
        }

        if (config.contains("eta_cut")){
                auto eta_cut= config["eta_cut"];
                eta_min = eta_cut[0];
                eta_max = eta_cut[1];
                if (eta_min <0.0 || eta_max <0.0 || eta_min>eta_max){
                        std::cout << "Invalid config file: interval 'eta_cut' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'eta_cut' key." << std::endl;
                return 1;
        }

        if (config.contains("nentries")){
                nentries = config["nentries"];
        }else{
                std::cerr << "Invalid config file: missing 'nentries' key." << std::endl;
                return 1;
        }

        if (config.contains("input_file")) {
                input_file = config["input_file"].template get<std::string>();
                if (!Config::exists0(input_file)){
                        std::cout << "Invalid config file: File " << input_file << " does not exist." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'input_file' key." << std::endl;
                return 1;
        }
        if (config.contains("DM_range")){
                auto DM_range = config["DM_range"];
                minDM = DM_range[0];
                maxDM = DM_range[1];
                if (minDM <0.0 || maxDM <0.0 || minDM>maxDM){
                        std::cout << "Invalid config file: interval 'DM_range' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'DM_range' key." << std::endl;
                return 1;
        }
        if (config.contains("BMcorr_range")){
                auto BMcorr_range = config["BMcorr_range"];
                minBMcorr = BMcorr_range[0];
                maxBMcorr = BMcorr_range[1];
                if (minBMcorr <0.0 || maxBMcorr <0.0 || minBMcorr>maxBMcorr){
                        std::cout << "Invalid config file: interval 'BMcorr_range' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'BMcorr_range' key." << std::endl;
                return 1;
        }
        if (config.contains("nvar_md")){
                nvar_md = config["nvar_md"];
                if (nvar_md<0){
                        std::cout << "Invalid config file: 'nvar_md' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nvar_md' key." << std::endl;
                return 1;
        }
         if (config.contains("nvar_mb")){
                nvar_mb = config["nvar_mb"];
                if (nvar_mb<0){
                        std::cout << "Invalid config file: 'nvar_mb' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nvar_mb' key." << std::endl;
                return 1;
        }
        if (config.contains("nbins")){
                nbins = config["nbins"];
                if (nbins<0){
                        std::cout << "Invalid config file: 'nbins' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nbins' key." << std::endl;
                return 1;
        }
	if (config.contains("DMshapes")) {
                DMshapes = config["DMshapes"].template get<std::vector<std::string>>();
                if (int(DMshapes.size())!=ncontr){
                        std::cout << "Invalid config file: Number of given DM pdfs does not match ncontr. " << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'DMshapes' key." << std::endl;
                return 1;
        }
	if (config.contains("BMshapes")) {
                BMshapes = config["BMshapes"].template get<std::vector<std::string>>();
                if (int(BMshapes.size())!=ncontr){
                        std::cout << "Invalid config file: Number of given BM pdfs does not match ncontr. " << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'BMshapes' key." << std::endl;
                return 1;
        }
	
	//mapping string values of shape vector onto integer values to use switch afterwards
        const std::unordered_map<std::string, const int> dictionaryDM ={
                {"DCBplusGaus", 0},
                {"Chebyshev", 1}
        };

        const std::unordered_map<std::string, const int> dictionaryBMcorr ={

                {"RCplusGaus", 0}
        };
	int iii=0;
	for (auto shape: DMshapes){
                if (dictionaryDM.find(shape) != dictionaryDM.end()) {
                        intshapesDM.push_back(dictionaryDM.at(shape));
			std::cout << shape << intshapesDM[iii++] << std::endl;

                } else {
                        std::cerr << "Shape '" << shape << "' not found." << std::endl;
                        return 1;
                }
        }
        std::cout << BMshapes.size() << std::endl;
        for (auto shape: BMshapes){
                if (dictionaryBMcorr.find(shape) != dictionaryBMcorr.end()) {
                        intshapesBMcorr.push_back(dictionaryBMcorr.at(shape));

                } else {
                        std::cerr << "Shape '" << shape << "' not found." << std::endl;
                        return 1;
                }
        }

	if (config.contains("fixVect")) {
                fixVect = config["fixVect"].template get<std::vector<std::string>>();
        } else {
                std::cerr << "Invalid config file: missing 'fixVect' key." << std::endl;
                return 1;
        }
	return 0;
}

