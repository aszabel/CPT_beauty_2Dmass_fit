#include "D_M_fit_shape.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis
{

	double DoubleSidedCrystalballPlusGaussPDF::EvalPDF(const double *xx, const double *par, const int component)
	{
		double f = abs(par[0]);
		double par1[2];
		double par2[6];
		par1[0] = par[1]; // mean
		par1[1] = par[2]; // sigma
		par2[0] = par[3]; // mean CB
		par2[1] = par[4]; // sigma CB
		par2[2] = par[5]; // alpha
		par2[3] = par[6]; // n
		par2[4] = par[7]; // alpha2
		par2[5] = par[8]; // n2

		if (component == 0)
			return (1 - f) * Gauss.EvalPDF(xx, par1);
		else if (component == 1)
			return f * DSCB.EvalPDF(xx, par2);
		else
			return (1 - f) * Gauss.EvalPDF(xx, par1) + f * DSCB.EvalPDF(xx, par2);
	}

	void DoubleSidedCrystalballPlusGaussPDF::CalcIntegral(const double *par, double min, double max)
	{
		double par1[2];
		double par2[6];
		par1[0] = par[1]; // mean
		par1[1] = par[2]; // sigma
		par2[0] = par[3]; // mean CB
		par2[1] = par[4]; // sigma CB
		par2[2] = par[5]; // alpha
		par2[3] = par[6]; // n
		par2[4] = par[7]; // alpha2
		par2[5] = par[8]; // n2

		Gauss.CalcIntegral(par1, min, max);
		DSCB.CalcIntegral(par2, min, max);
	}

	double JohnsonPlusGaussPDF::EvalPDF(const double *xx, const double *par, const int component)
	{
		double f = abs(par[0]);
		double par1[2];
		double par2[4];
		par1[0] = par[1]; // mean
		par1[1] = par[2]; // sigma
		par2[0] = par[3]; // xi
		par2[1] = par[4]; // lambda
		par2[2] = par[5]; // gamma
		par2[3] = par[6]; // delta

		if (component == 0)
			return (1 - f) * Gauss.EvalPDF(xx, par1);
		else if (component == 1)
			return f * JSU.EvalPDF(xx, par2);
		else
			return (1 - f) * Gauss.EvalPDF(xx, par1) + f * JSU.EvalPDF(xx, par2);
	}

	void JohnsonPlusGaussPDF::CalcIntegral(const double *par, double min, double max)
	{
		double par1[2];
		double par2[6];
		par1[0] = par[1]; // mean
		par1[1] = par[2]; // sigma
		par2[0] = par[3]; // xi
		par2[1] = par[4]; // lambda
		par2[2] = par[5]; // gamma
		par2[3] = par[6]; // delta

		Gauss.CalcIntegral(par1, min, max);
		JSU.CalcIntegral(par2, min, max);
	}

}
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
