#ifndef M_B_PDF_H
#define M_B_PDF_H
#include "pdf_interface.h"
#include <vector>

namespace cpt_b0_analysis{

class RaisedCosinePlusGaussPDF: public PDFInterface{
  	public:
	  	RaisedCosinePlusGaussPDF();
	  	~RaisedCosinePlusGaussPDF(){};

                void CalcIntegral(const double *par, double min, double max);
                double EvalPDF(const double *xx, const double *par);

	private:
		double IntCos;
		double IntGaus1;
		double IntGaus2;
};
}

#endif
