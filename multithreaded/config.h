// config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include "json.hpp"

class Config {
public:
	static int load(const std::string& filename);

	static inline bool exists0 (const std::string& name) {
		std::ifstream f(name.c_str());
    		return f.good();
    	}


    	static int sign;
    	static std::string input_file;

    	static int nentries;
    	static double event_fraction_for_correlation;
    	static double  event_fraction_final;

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

private:
    	Config() = default;
};

#endif // CONFIG_H
