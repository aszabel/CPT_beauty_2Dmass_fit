#ifndef KDE_FIT_H
#define KDE_FIT_H

namespace cpt_b0_analysis{

class mb_2misspt_fit{
  	public:
		const double minn, maxx; 
	  	mb_2misspt_fit(double _minn, double _maxx);
	  	~mb_2misspt_fit(){};

		double func_full(const double *xx, const double *par);
};
}

#endif
