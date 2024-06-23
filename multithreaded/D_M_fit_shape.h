#ifndef DM_FIT_H
#define DM_FIT_H
#include "pdf_interface.h"
#include <vector>

namespace cpt_b0_analysis{

class DoubleSidedCrystalballPlusGaussPDF: public PDFInterface{
  	public:
		
	  	DoubleSidedCrystalballPlusGaussPDF();
	  	~DoubleSidedCrystalballPlusGaussPDF(){};

		void CalcIntegral(const double *par, double min, double max);
		double EvalPDF(const double *xx, const double *par);

	private:
		double IntGaus;
		double IntDCB;
};
}

#endif
