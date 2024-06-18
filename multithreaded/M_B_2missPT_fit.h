#ifndef KDE_FIT_H
#define KDE_FIT_H
#include <vector>

namespace cpt_b0_analysis{

class mb_2misspt_fit{
  	public:
		const double minn, maxx; 
	  	mb_2misspt_fit(double _minn, double _maxx, int _ncontr);
	  	~mb_2misspt_fit(){};

		double func_full(const double *xx, const double *par, int icontr);
		int ncontr;
		std::vector<double> int_cos;
		std::vector<double> int_gaus1;
		std::vector<double> int_gaus2;
};
}

#endif
