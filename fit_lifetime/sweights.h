#ifndef SWEIGHTS_H
#define SWEIGHTS_H
#include <string>
#include <vector>

namespace cpt_b0_analysis
{
	class sWeights
	{
		public:
		sWeights(const std::string& _path_data);
		std::vector<std::pair<double, double>> get_sWeigths(double res[], bool BBbar);
		std::string path_data;

	};
		

}

#endif
