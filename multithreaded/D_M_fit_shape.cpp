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

	md_fit::md_fit(double _minx, double _maxx, int _ncontr) : minx(_minx), maxx(_maxx), ncontr(_ncontr)
	{
		for (int i = 0; i < ncontr; i++)
		{
			int_gaus.push_back(1.0);
			int_DCB.push_back(1.0);
		}
	}

	double md_fit::func_full(const double *xx, const double *_par, int icontr)
	{
		auto gausPDF = [this, icontr](const double *x, const double *par) -> double
		{
			double m_rec = x[0];
			double sigma = par[0];
			double mean = par[1];
			double f12 = abs(par[3]);
			double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
			if (int_gaus[icontr] != 0)
				gaus /= int_gaus[icontr];

			return f12 * gaus;
		};

		auto DoubleSidedCrystalballFunction = [this, icontr](const double *x, const double *par) -> double
		{
			double alpha = par[4];
			double n = par[5];
			double mean = par[1];
			double sigma = par[2];
			double f12 = abs(par[3]);
			double alpha_h = abs(par[6]);
			double m_rec = x[0];
			double result;

			if (m_rec < mean)
			{
				result = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
			}
			else
			{
				result = ROOT::Math::crystalball_function(2. * mean - m_rec, alpha_h, n, sigma, mean);
			}

			return (1.0 - f12) * result / int_DCB[icontr];
		};

		auto Deb_pdf = [this](const double *x, const double *par) -> double
		{
			double m_rec = x[0];
			double a1 = par[0];
			double a2 = par[1];
			double cheb = 1.0 + a1 * m_rec + a2 * (2.0 * m_rec * m_rec - 1.0);
			cheb /= int_Cheb;
			return cheb;
		};

		if (icontr != 4)
			return gausPDF(xx, _par) + DoubleSidedCrystalballFunction(xx, _par);
		else
			return Deb_pdf(xx, _par);
	}
}
