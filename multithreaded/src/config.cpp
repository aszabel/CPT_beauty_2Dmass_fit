// config.cpp
#include "config.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>

using json = nlohmann::json;
using namespace cpt_b0_analysis;


bool Config::isMC = false;
bool Config::binned = false;
bool Config::calc_sWeights = false;
bool Config::start_from_previous = false;
int Config::sign = -1;
std::string Config::input_file = "";
std::string Config::chainName = "";
std::string Config::previous_result_file = "";

int Config::nentries = -1;

std::vector<double> Config::tolerance = {};
int Config::functionCalls = 1e7;
int Config::printLevel = 1;
int Config::randSeed = -1;

double Config::muPTmin = 500.;
double Config::muPmin = 5000.;
double Config::eta_min  = 2.;
double Config::eta_max = 4.5;
double Config::tMin  = 1.0;
double Config::tMax = 15.0;

//Global
double Config::minDM = 1800.;
double Config::maxDM = 1940.;
double Config::minBMcorr = 2700.;
double Config::maxBMcorr = 8300.;
int Config::nvar_md = 7;
int Config::nvar_mb = 7;
int Config::ncontr = 6;
int Config::ntries = 1;

std::vector<int> Config::int_choose_fits = {};
std::vector<int> Config::intshapesDM = {};
std::vector<int> Config::intshapesBMcorr = {};

std::vector<std::string> Config::Fits = {};
std::vector<std::string> Config::DMshapes = {};
std::vector<std::string> Config::BMshapes = {};

std::vector<std::string> Config::fixVect = {};
std::vector<double> Config::fracInit = {};

std::vector<std::string> Config::contrName = {};
std::vector<std::string> Config::varname_md = {};
std::vector<std::string> Config::varname_mb = {};
std::map<std::string, std::string> Config::replace_var;
std::map<std::string, std::pair<double, double>> Config::varLimitsMap;

std::vector<std::vector<double>> Config::MC_MD={};
std::vector<std::vector<double>> Config::dMC_MD={};
std::vector<std::vector<double>> Config::MC_MB={};
std::vector<std::vector<double>> Config::dMC_MB = {};
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
void Config::read_MC(std::vector<std::vector<double>>& xx, std::vector<std::vector<double>>& dxx, std::string MC_directory, int nvar){
// Read results of 1D fits that are stored as simple text files
        // Each line corresponds to a single parameter
        // First column defines parameter value
        // Second column defines parameter uncertainty
        // D mass fit
	xx.reserve(ncontr);
	dxx.reserve(ncontr);
        for (int ifile = 0; ifile < ncontr; ifile++)
        {
		const int nbuffer = 100;
		char name[nbuffer];
		snprintf(name, nbuffer, "%s/res_%s_%d.txt", MC_directory.c_str(), contrName[ifile].c_str(), sign);
                std::ifstream infile(name);
		double _x, _dx;
		std::vector<double> tmp_x = {};
		std::vector<double> tmp_dx = {};
                while (infile >> _x >> _dx ){
			tmp_x.push_back(_x);
			tmp_dx.push_back(_dx);
		}

  		xx.emplace_back(std::move(tmp_x));
  		dxx.emplace_back(std::move(tmp_dx));
		if (int(xx[ifile].size())!=nvar || int(dxx[ifile].size())!=nvar){
			std::cerr << " Wrong number of parameters in the MC result file " << name << std::endl;
        	}
	}
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
                        std::cerr << "Invalid config file: allowed values for the 'sign' key are 'muplus' or 'muminus'." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'sign' key." << std::endl;
                return 1;
        }

        if (config.contains("isMC")) {
                if (config["isMC"].template get<std::string>() == std::string("true"))
                        isMC = true;
                else if (config["isMC"].template get<std::string>() == std::string("false"))
                        isMC = false;
                else {
                        std::cerr << "Invalid config file: allowed values for the 'isMC' key are 'true' or 'false'." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'isMC' key." << std::endl;
                return 1;
        }

        if (config.contains("binned")) {
                if (config["binned"].template get<std::string>() == std::string("true"))
                        binned = true;
                else if (config["binned"].template get<std::string>() == std::string("false"))
                        binned = false;
                else {
                        std::cerr << "Invalid config file: allowed values for the 'binned' key are 'true' or 'false'." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'binned' key." << std::endl;
                return 1;
        }

        if (config.contains("calc_sWeights")) {
                if (config["calc_sWeights"].template get<std::string>() == std::string("true"))
                        calc_sWeights = true;
                else if (config["calc_sWeights"].template get<std::string>() == std::string("false"))
                        calc_sWeights = false;
                else {
                        std::cerr << "Invalid config file: allowed values for the 'calc_sWeights' key are 'true' or 'false'." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'calc_sWeights' key." << std::endl;
                return 1;
        }

        if (config.contains("start_from_previous")) {
                if (config["start_from_previous"].template get<std::string>() == std::string("true"))
                        start_from_previous = true;
                else if (config["start_from_previous"].template get<std::string>() == std::string("false"))
                        start_from_previous = false;
                else {
                        std::cerr << "Invalid config file: allowed values for the 'start_from_previous' key are 'true' or 'false'." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'start_from_previous' key." << std::endl;
                return 1;
        }

	struct stat sb;
       	if (start_from_previous) {
                if (config.contains("previous_result_file")) {
                        previous_result_file = config["previous_result_file"].template get<std::string>();
                        if (stat(previous_result_file.c_str(), &sb)!=0 && start_from_previous){
                                std::cerr << "Invalid config file: Previous file " <<  previous_result_file << " does not exist." << std::endl;
                                return 1;
                        }
                } else {
                        std::cerr << "Invalid config file: missing 'previous_result_file' key." << std::endl;
                        return 1;
                }
        }
 

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


        if (config.contains("nentries")){
                nentries = config["nentries"];
        }else{
                std::cerr << "Invalid config file: missing 'nentries' key." << std::endl;
                return 1;
        }

        if (config.contains("input_file")) {
                input_file = config["input_file"].template get<std::string>();
                if (stat(input_file.c_str(), &sb)!=0){
                        std::cerr << "Invalid config file: Input file " << input_file << " does not exist." << std::endl;
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
	
	if (config.contains("Fits")) {
                Fits = config["Fits"].template get<std::vector<std::string>>();
                if (int(Fits.size())==ncontr){
                        std::cerr << "Invalid config file: no fits chosen " << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'Fits' key." << std::endl;
                return 1;
        }	 
	for (auto fit: Fits){
                if (dictionaryChooseFit.find(fit) != dictionaryChooseFit.end()) {
                        int_choose_fits.push_back(dictionaryChooseFit.at(fit));

                } else {
                        std::cerr << "Fit '" << fit << "' not found. Chose 'frac', 'BM', 'DM+BMfixed', 'shapes' or 'all'" << std::endl;
                        return 1;
                }
        }
	if (config.contains("tolerance")) {
                tolerance = config["tolerance"].template get<std::vector<double>>();
                if(tolerance.size()!=Fits.size()){
                        std::cerr << " You have to provide the same number of Fits and tolerances\n";
                        return 1;
                }
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
                fixVect = config["fixVect"].template get<std::vector<std::string>>();
        } else {
                std::cerr << "Invalid config file: missing 'fixVect' key." << std::endl;
                return 1;
        }

	
	if (config.contains("fracInit")) {
                fracInit = config["fracInit"].template get<std::vector<double>>();
		if( int(fracInit.size())!=ncontr){
			std::cerr << "Invalid config file: vector 'fracInit' should have " << ncontr << " elements.";
		        return 1;	
		}
        } else {
                std::cerr << "Invalid config file: missing 'fracInit' key." << std::endl;
                return 1;
        }
        if (config.contains("contrName")) {
                contrName = config["contrName"].template get<std::vector<std::string>>();
		if( int(contrName.size())!=ncontr){
			std::cerr << "Invalid config file: vector 'contrName' should have " << ncontr << " elements.";
		        return 1;	
		}
        } else {
                std::cerr << "Invalid config file: missing 'contrName' key." << std::endl;
                return 1;
        }
	if (config.contains("varname_md")) {
                varname_md = config["varname_md"].template get<std::vector<std::string>>();
                if( int(varname_md.size())!=nvar_md){
                        std::cerr << "Invalid config file: vector 'varname_md' should have " << nvar_md << " elements.";
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'varname_md' key." << std::endl;
                return 1;
        }
	if (config.contains("varname_mb")) {
                varname_mb = config["varname_mb"].template get<std::vector<std::string>>();
                if( int(varname_mb.size())!=nvar_mb){
                        std::cerr << "Invalid config file: vector 'varname_mb' should have " << nvar_mb<< " elements.";
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'varname_mb' key." << std::endl;
                return 1;
        }

	std::string MC_directory_MD;
        if (config.contains("MC_directory_MD")) {
                MC_directory_MD = config["MC_directory_MD"].template get<std::string>();
                if (stat(MC_directory_MD.c_str(), &sb)!=0){
                        std::cerr << "Invalid config file: MC_directory_MD " << MC_directory_MD << " does not exist." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'MC_directory_MD' key." << std::endl;
                return 1;
        }
	std::string MC_directory_MB;
        if (config.contains("MC_directory_MB")) {
                MC_directory_MB= config["MC_directory_MB"].template get<std::string>();
                if (stat(MC_directory_MB.c_str(), &sb)!=0){
                        std::cerr << "Invalid config file: MC_directory_MB " << MC_directory_MB << " does not exist." << std::endl;
                        return 1;
                }
        } else {
                std::cerr << "Invalid config file: missing 'MC_directory_MB' key." << std::endl;
                return 1;
        }

        if(!isMC) {
                // Initial fit parameter values taken from 1D fits to MC and Side Bands
                read_MC(MC_MD, dMC_MD, MC_directory_MD, nvar_md);
                read_MC(MC_MB, dMC_MB, MC_directory_MB, nvar_mb);
        }

        if (config.contains("replace_var")) {
                json entry = config["replace_var"];
                for (auto it = entry.begin(); it!=entry.end();++it){
                        std::string aim = it.value().template get<std::string>();
                        replace_var[it.key()] = aim;
		}
                for (const auto& rep_var: replace_var){
                        fixVect.push_back(rep_var.first);
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

