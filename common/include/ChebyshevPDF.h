#ifndef CHEB_PDF_H
#define CHEB_PDF_H
#include "pdf_interface.h"

namespace cpt_b0_analysis
{

	class ChebyshevPDF : public PDFInterface
	{
	public:
		ChebyshevPDF();
		ChebyshevPDF(const ChebyshevPDF&) = default;
		ChebyshevPDF(ChebyshevPDF&&) = default;
		~ChebyshevPDF() {};

		void CalcIntegral(const double *par, double min, double max);
		double EvalPDF(const double *xx, const double *par, const int component=-1);
		double getIntegral();
		int getComponentCount(){return 0;};

	private:
		double IntCheb;
	};


	class ExponentPDF : public PDFInterface
	{
	public:
		ExponentPDF();
		ExponentPDF(const ExponentPDF&) = default;
		ExponentPDF(ExponentPDF&&) = default;
		~ExponentPDF() {};

		void CalcIntegral(const double *par, double min, double max);
		double EvalPDF(const double *xx, const double *par, const int component=-1);
		double getIntegral();
		int getComponentCount(){return 0;};

	private:
		double IntExp;
	};
}

#endif
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
