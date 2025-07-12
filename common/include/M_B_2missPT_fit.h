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
		double EvalPDF(const double *xx, const double *par);

	private:
		double IntCos;
		double IntGaus1;
		double IntGaus2;
	};

	class SkewNormalPlusGausPDF : public PDFInterface
	{
	public:
		SkewNormalPlusGausPDF();
		SkewNormalPlusGausPDF(const SkewNormalPlusGausPDF &other)=default;
		SkewNormalPlusGausPDF(SkewNormalPlusGausPDF&&)=default;
		~SkewNormalPlusGausPDF() {};
			
		void CalcIntegral(const double *par, double min, double max);
		double getIntegral(){return 0.0;};
		double EvalPDF(const double *xx, const double *par);

	private:
		double IntSkewNorm;
		double IntCB;
		std::function<double(const double*, const double*)> skew_normal;
	};
}
	

#endif
