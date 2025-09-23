#include "M_fit_Every.h"

using namespace cpt_b0_analysis;
void B_M_fit_Every(std::string config_file){
	// Load config
	std::cout<<config_file<<endl;
	if (Config::load(config_file)){
		std::cerr<< " Bad config file! " << std::endl;
		return;
	}

	TString binning = "";
	if (Config::binned)
		binning = "binned";
	else
		binning = "unbinned";

	TString path_results = Form("%s_%s", Config::MC_directory_MB.c_str(), binning.Data());

	M_fit_Every(
		path_results,
		Config::nvar_mb,
		Config::minBMcorr,
		Config::maxBMcorr,
		TString("Bmass"),
		Config::varname_mb
	);
}

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
