#ifndef CHEB_PDF_H
#define CHEB_PDF_H
#include "pdf_interface.h"

namespace cpt_b0_analysis
{

	class ChebyshevPDF : public PDFInterface
	{
	public:
		ChebyshevPDF();
		~ChebyshevPDF() {};

		void CalcIntegral(const double *par, double min, double max);
		double EvalPDF(const double *xx, const double *par);

	private:
		double IntCheb;
	};
}

#endif
