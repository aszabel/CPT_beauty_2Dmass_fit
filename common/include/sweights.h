#ifndef SWEIGHTS_H
#define SWEIGHTS_H
#include <string>
#include <vector>

namespace cpt_b0_analysis
{
	class sWeights
	{
		public:
		sWeights(const std::string& _path_data, const std::string& _outfilename, const std::string& _outtreename);
		void get_sWeigths(const double res[], bool BBbar);
		std::string path_data, outfilename, outtreename;

	};
		

}

#endif
