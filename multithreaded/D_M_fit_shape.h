#ifndef DM_FIT_H
#define DM_FIT_H
#include <vector>

namespace cpt_b0_analysis
{

	class md_fit
	{
	public:
		const double minx, maxx;
		md_fit(double _minn, double _maxx, int _ncontr);
		~md_fit(){};

		double func_full(const double *xx, const double *par, int icontr);
		int nevents;
		int ncontr;
		std::vector<double> int_gaus;
		std::vector<double> int_DCB;
		double int_Cheb;
	};
}

#endif
