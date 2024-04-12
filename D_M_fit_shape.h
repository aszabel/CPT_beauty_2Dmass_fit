#ifndef DM_FIT_H
#define DM_FIT_H

namespace cpt_b0_analysis{

class md_fit{
  	public:
		const double minx, maxx; 
	  	md_fit(double _minn, double _maxx);
	  	~md_fit(){};

		double func_full(const double *xx, const double *par);
		int nevents;
};
}

#endif
