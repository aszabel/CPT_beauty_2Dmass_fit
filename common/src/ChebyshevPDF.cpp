#include "ChebyshevPDF.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis
{

	ChebyshevPDF::ChebyshevPDF()
	{
		IntCheb = 1.0;
	}
	double ChebyshevPDF::EvalPDF(const double *xx, const double *_par, const int component)
	{
		auto DebPDF = [this](const double *x, const double *par) -> double
		{
			double m_rec = x[0];
			double a1 = par[0];
			double a2 = par[1];
			double cheb = 1.0 + a1 * m_rec + a2 * (2.0 * m_rec * m_rec - 1.0);
			cheb /= IntCheb;
			return cheb;
		};

		return DebPDF(xx, _par);
	}
	void ChebyshevPDF::CalcIntegral(const double *par, double min, double max)
	{
		double a1 = par[0];
		double a2 = par[1];

		IntCheb = (1.0 - a2) * max + 0.5 * a1 * max * max + 2. / 3. * a2 * max * max * max - (1.0 - a2) * min - 0.5 * a1 * min * min - 2. / 3. * a2 * min * min * min;
	}
	double ChebyshevPDF::getIntegral(){
		return IntCheb;
	}

	ExponentPDF::ExponentPDF()
	{
		IntExp = 1.0;
	}
	double ExponentPDF::EvalPDF(const double *xx, const double *_par, const int component)
	{
		auto ExpPDF = [this](const double *x, const double *par) -> double
		{
			double m_rec = x[0];
			double slope = par[0];
			double exp = TMath::Exp(-slope*m_rec);
			if (IntExp!=0.0)
				exp /= IntExp;
			return exp;
		};

		return ExpPDF(xx, _par);
	}
	void ExponentPDF::CalcIntegral( const double *par, double min, double max)
	{
		double slope = par[0];
		if(slope!=0.0)
			IntExp = -1./slope*TMath::Exp(-slope*max)+1./slope*TMath::Exp(-slope*min);
		else 
			IntExp = 1.0;
	}

	double ExponentPDF::getIntegral(){
		return IntExp;
	}
}
