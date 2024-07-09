#ifndef PDF_INTERFACE_H
#define PDF_INTERFACE_H

namespace cpt_b0_analysis
{

	class PDFInterface
	{
	public:
		PDFInterface() {}
		PDFInterface(const PDFInterface& other) {}
		virtual ~PDFInterface() {}
		virtual void CalcIntegral(const double *par, double min, double max) = 0; // "= 0" part makes this method pure virtual, and also makes this class abstract.
		virtual double EvalPDF(const double *par, const double *x) = 0;
	};

}
#endif
