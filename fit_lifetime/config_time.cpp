// config.cpp
#include "config_time.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>

using json = nlohmann::json;
using namespace cpt_b0_analysis;


std::string Config::input_file = "";
std::string Config::r_factor_file = "";
std::string Config::chainName = "";

int Config::nentries = -1;

double Config::tolerance = 0.1;
int Config::functionCalls = 1e7;
int Config::printLevel = 1;
int Config::randSeed = -1;

double Config::muPTmin = 500.;
double Config::muPmin = 5000.;
double Config::eta_min  = 2.;
double Config::eta_max = 4.5;
double Config::tMin  = 1.0;
double Config::tMax = 15.0;
double Config::r_factor_min  = 0.0;
double Config::r_factor_max = 2.0;
double Config::draw_min  = 0.99;
double Config::draw_max = 1.01;


double Config::Rez = 0.0;
double Config::Imz = 0.0;
double Config::eta = 0.99728;
double Config::alpha = 9.97901;
double Config::beta = -0.13718;
double Config::t_shift = -1.24772;



double Config::imz_start = 0.0;

//Global
double Config::minDM = 1800.;
double Config::maxDM = 1940.;
double Config::minBMcorr = 2700.;
double Config::maxBMcorr = 8300.;
int Config::Nbins = 250;
int Config::nlogBins = 50;
int Config::nvar_md = 7;
int Config::nvar_mb = 7;
int Config::ncontr = 6;
int Config::ntries = 1;

std::vector<int> Config::intshapesDM = {};
std::vector<int> Config::intshapesBMcorr = {};

std::vector<std::string> Config::DMshapes = {};
std::vector<std::string> Config::BMshapes = {};


std::vector<std::vector<std::string>> Config::fixVect = {};

std::vector<std::pair<std::string, double>> Config::varname;
int Config::nvar_time = 10;
std::map<std::string, std::string> Config::replace_var;
std::map<std::string, std::pair<double, double>> Config::varLimitsMap;

std::vector<std::shared_ptr<PDFInterface>> Config::getVectorPDFs(const std::string& domain) {
	std::vector<std::shared_ptr<PDFInterface>> vPDFs;
	if (domain == std::string("Dmass")){
		for (auto intshape: intshapesDM){
                	switch (intshape){
                        	case 0:
                                	vPDFs.push_back(std::make_shared<DoubleSidedCrystalballPlusGaussPDF>());
                                	break;
                        	case 1:
                                	vPDFs.push_back(std::make_shared<ChebyshevPDF>());
                                	break;
                        	default:
                                	std::cerr << "Error while DM_pdf dynamic declaration. Check if shapes from config file refer to the shapes defined in the code." << std::endl;
                                	return std::vector<std::shared_ptr<PDFInterface>>{};
                	}
		}

	}else if (domain == std::string("Bmass")){
	        for (auto intshape: intshapesBMcorr){
        	        switch (intshape){
                        	case 0:
                                	vPDFs.push_back(std::make_shared<RaisedCosinePlusGaussPDF>());
                                	break;
                        	default:
                                	std::cerr << "Error while BM_pdf dynamic declaration. Check if shapes from config file refer to the shapes defined in the code." << std::endl;
                                return std::vector<std::shared_ptr<PDFInterface>>{};
                	}
        	}

	}else{
		std::cerr << "Error wrong domaine in getVectorPDFs(domain) chose 'Dmass' or 'Bmass'" << std::endl;
		return std::vector<std::shared_ptr<PDFInterface>>{};
	}	
    	return vPDFs;
}

int Config::load(const std::string& filename) {
	std::ifstream config_file(filename);
        json config = json::parse(config_file);
        config_file.close();


	if (config.contains("chainName")) {
                chainName  = config["chainName"].template get<std::string>();
        } else {
                std::cerr << "Invalid config file: missing 'chainName' key." << std::endl;
                return 1;
        }

        if (config.contains("functionCalls")){
                functionCalls = config["functionCalls"];
                if (functionCalls <0.){
                        std::cerr << "Invalid config file: 'functionCalls' must be a positive integer." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'functionCalls' key." << std::endl;
                return 1;
        }


        if (config.contains("printLevel")){
                printLevel = config["printLevel"];
                if (printLevel<0 || printLevel>3){
                        std::cerr << "Invalid config file: 'printLevel' must be 0, 1, 2 or 3" << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'printLevel' key." << std::endl;
                return 1;
        }

	if (config.contains("randSeed")){
                randSeed = config["randSeed"];
                if (randSeed<-1){
                        std::cerr << "Invalid config file: 'randSeed' must be -1, 0 or a positive integer." << std::endl; //when -1 then no randomization
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'randSeed' key." << std::endl;
                return 1;
        }


        if (config.contains("muPmin")){
                muPmin = config["muPmin"];
                if (muPmin <0){
                        std::cerr << "Invalid config file: 'muPmin' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'muPmin' key." << std::endl;
                return 1;
        }

        if (config.contains("muPTmin")){
                muPTmin = config["muPTmin"];
                if (muPmin <0){
                        std::cerr << "Invalid config file: 'muPTmin' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'muPTmin' key." << std::endl;
                return 1;
        }


        if (config.contains("eta_cut")){
                auto eta_cut= config["eta_cut"];
                eta_min = eta_cut[0];
                eta_max = eta_cut[1];
                if (eta_min <0.0 || eta_max <0.0 || eta_min>eta_max){
                        std::cerr << "Invalid config file: interval 'eta_cut' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'eta_cut' key." << std::endl;
                return 1;
        }

        if (config.contains("draw_cut")){
                auto draw_cut= config["draw_cut"];
                draw_min = draw_cut[0];
                draw_max = draw_cut[1];
        }else{
                std::cerr << "Invalid config file: missing 'draw_cut' key." << std::endl;
                return 1;
        }

        if (config.contains("t_cut")){
                auto t_cut= config["t_cut"];
                tMin = t_cut[0];
                tMax = t_cut[1];
                if (tMin <0.0 || tMax <0.0 || tMin>tMax){
                        std::cerr << "Invalid config file: interval 't_cut' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 't_cut' key." << std::endl;
                return 1;
        }

	if (config.contains("imz_start")){
               	imz_start = config["imz_start"];
        }else{
                std::cerr << "Invalid config file: missing 'imz_start' key." << std::endl;
                return 1;
        }


        if (config.contains("nentries")){
                nentries = config["nentries"];
        }else{
                std::cerr << "Invalid config file: missing 'nentries' key." << std::endl;
                return 1;
        }

	struct stat sb;
        if (config.contains("input_file")) {
                input_file = config["input_file"].template get<std::string>();
                if (stat(input_file.c_str(), &sb)!=0){
                        std::cerr << "Invalid config file: File " << input_file << " does not exist." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'input_file' key." << std::endl;
                return 1;
        }
       
        if (config.contains("r_factor_file")) {
                r_factor_file = config["r_factor_file"].template get<std::string>();
                if (stat(r_factor_file.c_str(), &sb)!=0){
                        std::cerr << "Invalid config file: File " << r_factor_file << " does not exist." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'r_factor_file' key." << std::endl;
                return 1;
        }
	
       if (config.contains("DM_range")){
                auto DM_range = config["DM_range"];
                minDM = DM_range[0];
                maxDM = DM_range[1];
                if (minDM <0.0 || maxDM <0.0 || minDM>maxDM){
                        std::cerr << "Invalid config file: interval 'DM_range' must be positive." << std::endl;
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
                        std::cerr << "Invalid config file: interval 'BMcorr_range' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'BMcorr_range' key." << std::endl;
                return 1;
        }
        if (config.contains("r_factor_range")){
                auto r_factor_range = config["r_factor_range"];
                r_factor_min = r_factor_range[0];
                r_factor_max = r_factor_range[1];
        }else{
                std::cerr << "Invalid config file: missing 'r_factor_range' key." << std::endl;
                return 1;
        }
        if (config.contains("Nbins")){
                Nbins = config["Nbins"];
                if (Nbins<0){
                        std::cerr << "Invalid config file: 'Nbins' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'Nbins' key." << std::endl;
                return 1;
        }        
	
	if (config.contains("nlogBins")){
                nlogBins = config["nlogBins"];
                if (nlogBins<0){
                        std::cerr << "Invalid config file: 'nlogBins' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nlogBins' key." << std::endl;
                return 1;
        }
 
         if (config.contains("nvar_md")){
                nvar_md = config["nvar_md"];
                if (nvar_md<0){
                        std::cerr << "Invalid config file: 'nvar_md' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nvar_md' key." << std::endl;
                return 1;
        }
        if (config.contains("nvar_mb")){
                nvar_mb = config["nvar_mb"];
                if (nvar_mb<0){
                        std::cerr << "Invalid config file: 'nvar_mb' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nvar_mb' key." << std::endl;
                return 1;
        } if (config.contains("nvar_time")){
                nvar_time = config["nvar_time"];
                if (nvar_mb<0){
                        std::cerr << "Invalid config file: 'nvar_time' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'nvar_time' key." << std::endl;
                return 1;
        }
	         if (config.contains("ntries")){
                ntries = config["ntries"];
                if (ntries<=0){
                        std::cerr << "Invalid config file: 'ntries' must be positive." << std::endl;
                        return 1;
                }
        }else{
                std::cerr << "Invalid config file: missing 'ntries' key." << std::endl;
                return 1;
        }
	
	if (config.contains("tolerance")) {
                tolerance = config["tolerance"].template get<double>();
        } else {
                std::cerr << "Invalid config file: missing 'tolerance' key." << std::endl;
                return 1;
        }

	if (config.contains("DMshapes")) {
                DMshapes = config["DMshapes"].template get<std::vector<std::string>>();
                if (int(DMshapes.size())!=ncontr){
                        std::cerr << "Invalid config file: Number of given DM pdfs does not match ncontr. " << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'DMshapes' key." << std::endl;
                return 1;
        }
	if (config.contains("BMshapes")) {
                BMshapes = config["BMshapes"].template get<std::vector<std::string>>();
                if (int(BMshapes.size())!=ncontr){
                        std::cerr << "Invalid config file: Number of given BM pdfs does not match ncontr. " << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'BMshapes' key." << std::endl;
                return 1;
        }

	for (auto shape: DMshapes){
                if (dictionaryDM.find(shape) != dictionaryDM.end()) {
                        intshapesDM.push_back(dictionaryDM.at(shape));

                } else {
                        std::cerr << "Shape '" << shape << "' not found." << std::endl;
                        return 1;
                }
        }
        for (auto shape: BMshapes){
                if (dictionaryBMcorr.find(shape) != dictionaryBMcorr.end()) {
                        intshapesBMcorr.push_back(dictionaryBMcorr.at(shape));

                } else {
                        std::cerr << "Shape '" << shape << "' not found." << std::endl;
                        return 1;
                }
        }

	if (config.contains("fixVect")) {
                fixVect = config["fixVect"].template get<std::vector<std::vector<std::string>>>();
        } else {
                std::cerr << "Invalid config file: missing 'fixVect' key." << std::endl;
                return 1;
        }

	
        if (config.contains("varname")) {
                varname = config["varname"].template get<std::vector<std::pair<std::string, double>>>();
                //for (auto it = entry.begin(); it!=entry.end();++it){
                        //double target = it.value().template get<double>();
                        //varname[it.key()] = target;
			//varname.insert(*it);
		//}
        } else {
                std::cerr << "Invalid config file: missing 'varname' key." << std::endl;
                return 1;
        }

        if (config.contains("replace_var")) {
                json entry = config["replace_var"];
                for (auto it = entry.begin(); it!=entry.end();++it){
                        std::string target = it.value().template get<std::string>();
                        replace_var[it.key()] = target;
		}
        } else {
                std::cerr << "Invalid config file: missing 'replace_var' key." << std::endl;
                return 1;
        }

        if (config.contains("varLimitsVect")) {
                json entry = config["varLimitsVect"];
                for( auto it = entry.begin(); it!=entry.end();++it){
                        std::pair<double, double> pair_lims = it.value().template get<std::pair<double, double>>();
                        varLimitsMap[it.key()] = pair_lims;
                }
        } else {
                std::cerr << "Invalid config file: missing 'varLimitsVect' key." << std::endl;
                return 1;
        }


	return 0;
}

