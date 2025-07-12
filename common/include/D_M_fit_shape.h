#ifndef DM_FIT_H
#define DM_FIT_H
#include "pdf_interface.h"
#include <vector>

namespace cpt_b0_analysis
{

	class DoubleSidedCrystalballPlusGaussPDF : public PDFInterface
	{
	public:
		DoubleSidedCrystalballPlusGaussPDF();
		DoubleSidedCrystalballPlusGaussPDF(const DoubleSidedCrystalballPlusGaussPDF& other)=default;
		DoubleSidedCrystalballPlusGaussPDF(DoubleSidedCrystalballPlusGaussPDF&&)=default;
		~DoubleSidedCrystalballPlusGaussPDF() {};

		void CalcIntegral(const double *par, double min, double max);
		double getIntegral(){return 0.0;};
		double EvalPDF(const double *xx, const double *par);

	private:
		double IntGaus;
		double IntDCB;
	};


	class JohnsonPlusGaussPDF : public PDFInterface
	{
	public:
		JohnsonPlusGaussPDF();
		JohnsonPlusGaussPDF(const JohnsonPlusGaussPDF& other)=default;
		JohnsonPlusGaussPDF(JohnsonPlusGaussPDF&&)=default;
		~JohnsonPlusGaussPDF() {};

		void CalcIntegral(const double *par, double min, double max);
		double getIntegral(){return 0.0;};
		double EvalPDF(const double *xx, const double *par);

	private:
		double IntGaus;
		double IntJSU;
	};
}

#endif
