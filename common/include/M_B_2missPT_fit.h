#ifndef M_B_PDF_H
#define M_B_PDF_H
#include "pdf_interface.h"
#include <vector>
#include <functional>

namespace cpt_b0_analysis
{

	class RaisedCosinePlusGaussPDF : public PDFInterface
	{
	public:
		RaisedCosinePlusGaussPDF();
		RaisedCosinePlusGaussPDF(const RaisedCosinePlusGaussPDF &other)=default;
		RaisedCosinePlusGaussPDF(RaisedCosinePlusGaussPDF&&)=default;
		~RaisedCosinePlusGaussPDF() {};

		void CalcIntegral(const double *par, double min, double max);
		double getIntegral(){return 0.0;};
		int getComponentCount(){return 3;}; // RaisedCosine + 2 x Gauss
		double EvalPDF(const double *xx, const double *par, const int component=-1);

	private:
		double IntCos;
		double IntGaus1;
		double IntGaus2;
	};

	class SkewNormalPlusCBPDF : public PDFInterface
	{
	public:
		SkewNormalPlusCBPDF();
		SkewNormalPlusCBPDF(const SkewNormalPlusCBPDF &other)=default;
		SkewNormalPlusCBPDF(SkewNormalPlusCBPDF&&)=default;
		~SkewNormalPlusCBPDF() {};
			
		void CalcIntegral(const double *par, double min, double max);
		double getIntegral(){return 0.0;};
		int getComponentCount(){return 2;}; // SkewNormal + CB
		double EvalPDF(const double *xx, const double *par, const int component = -1);

	private:
		double IntSkewNorm;
		double IntCB;
		double skew_normal(const double*, const double*);
	};
}

#endif
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
