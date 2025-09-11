#include "M_B_2missPT_fit.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include <Math/DistFunc.h>  // normal_cdf
#include <boost/math/special_functions/owens_t.hpp>

namespace cpt_b0_analysis
{

	RaisedCosinePlusGaussPDF::RaisedCosinePlusGaussPDF()
	{
		IntCos = 1.0;
		IntGaus1 = 1.0;
		IntGaus2 = 1.0;
	}

	double RaisedCosinePlusGaussPDF::EvalPDF(const double *xx, const double *par, const int component)
	{
		auto raised_cosine = [this](const double *x, const double *par) -> double
		{
			double mass = x[0];
			double s = abs(par[1]);
			double mean = par[0];
			double cos = 0.0;
			if (mass >= mean - s && mass <= mean + s)
			{
				cos = 1.0 / (2.0 * s) * (1.0 + TMath::Cos((mass - mean) / s * TMath::Pi()));
			}
			if (IntCos != 0.0)
				cos /= IntCos;
			return cos;
		};
		auto gaus = [this](const double *x, const double *par, const int idx) -> double
		{
			const double m_rec = x[0];
			double norm = 0.0;
			if (idx == 0)
				norm = 1.0 - abs(par[2]);
			else
				norm = abs(par[2]);
			if (norm<=1.0e-9) 
				norm = 0.0;

			double mean = abs(par[3 + idx * 2]);
			double sigma = abs(par[4 + idx * 2]);
			double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
			gaus *= norm;

			if (idx == 0) {
				if (abs(IntGaus1) >= 1.0e-6)
					gaus /= IntGaus1;
			} else {
				if (IntGaus2 != 0.0)
					gaus /= IntGaus2;
			}

			return gaus;
		};
		double frac = par[7];

		if (component == 0)
			return (1.0-frac)*raised_cosine(xx, par);
		else if (component == 1)
			return frac*gaus(xx, par, 0);
		else if (component == 2)
			return frac*gaus(xx, par, 1);
		else
			return ((1.0-frac)*raised_cosine(xx, par) + frac*(gaus(xx, par, 0) + gaus(xx, par, 1)));
	}

	void RaisedCosinePlusGaussPDF::CalcIntegral(const double *par, double min, double max)
	{
		double s = abs(par[1]);
		double mean = par[0];
		double minnB = mean - s;
		if (min > mean - s)
			minnB = min;
		double maxxB = mean + s;
		if (max < mean + s)
			maxxB = max;

		IntCos = 1.0 / (2.0) * (1.0 + (maxxB - mean) / s + TMath::Sin((maxxB - mean) / s * TMath::Pi()) / TMath::Pi()) - 1.0 / (2.0) * (1.0 + (minnB - mean) / s + TMath::Sin((minnB - mean) / s * TMath::Pi()) / TMath::Pi());

		double mean1 = abs(par[3]);
		double sigma1 = abs(par[4]);

		IntGaus1 = ROOT::Math::normal_cdf(max, sigma1, mean1) - ROOT::Math::normal_cdf(min, sigma1, mean1);

		if (abs(par[2])-1.0 < 1.0e-6){
			IntGaus1 = 0.0;
		} 
		double mean2 = abs(par[5]);
		double sigma2 = abs(par[6]);

		IntGaus2 = ROOT::Math::normal_cdf(max, sigma2, mean2) - ROOT::Math::normal_cdf(min, sigma2, mean2);
		if (abs(par[2])<1.0e-6) 
			IntGaus2 = 0.0;
	}

	SkewNormalPlusCBPDF::SkewNormalPlusCBPDF()
	{
		IntSkewNorm = 1.0;
		IntCB = 1.0;
	}
	
	double SkewNormalPlusCBPDF::skew_normal(const double *x, const double *par)
	{
		double xm = x[0]-par[6];
		double sigma = par[4];
		double skew  = par[5];

		if (skew == 0.0) return ROOT::Math::gaussian_pdf(xm, sigma, 0.0);

		return 2.*ROOT::Math::gaussian_pdf(xm, sigma, 0.0)*ROOT::Math::normal_cdf(skew*xm/sigma,1.0, 0.0);
	};

	double SkewNormalPlusCBPDF::EvalPDF(const double *xx, const double *par, const int component)
	{
		auto crystal_ball = [this](const double *x, const double *par) -> double
		{
			double m_rec = x[0];
			double sigma = par[0];
			double alpha = abs(par[2])+1.0e-6;
			double mean = par[1];
			double n = par[3];
			double CB = ROOT::Math::crystalball_function(2.0*mean-m_rec, alpha, n, sigma, mean);
			if (IntCB != 0.0)
				CB /= IntCB;
			return CB;
		};
		if (IntSkewNorm == 0.0)
			return 0.0;

		if (component == 0)
			return abs(1.0-par[7])*skew_normal(xx, par)/IntSkewNorm;
		else if (component == 1)
			return abs(par[7])*crystal_ball(xx, par);
		else
			return abs(1.0-par[7])*skew_normal(xx, par)/IntSkewNorm+abs(par[7])*crystal_ball(xx, par);
	}

	void SkewNormalPlusCBPDF::CalcIntegral(const double *par, double min, double max)
	{
		double xmin = min-par[6];
		double xmax = max-par[6];
		double sigma_sk = par[4];
		double skew  = par[5];

		double zmin = xmin / sigma_sk;
		double zmax = xmax / sigma_sk;
		double Phi = ROOT::Math::normal_cdf(zmax)-ROOT::Math::normal_cdf(zmin);

		if (skew == 0) {
			IntSkewNorm = Phi;
		} else {
			double T = boost::math::owens_t(zmax, skew)-boost::math::owens_t(zmin, skew);
			IntSkewNorm = Phi - 2.0 * T;
		}

/*
		TF1 skewfunc("skewfunc", skew_normal, min, max, 8);
		skewfunc.SetParameters(par);
		ROOT::Math::WrappedTF1 wf1(skewfunc);
		ROOT::Math::GaussIntegrator ig;
		ig.SetFunction(wf1);
		ig.SetRelTolerance(0.01);
		IntSkewNorm = skewfunc.Integral(min, max, 1.0e-8);
		if (!std::isfinite(IntSkewNorm)){
			for (int i=0; i<8; i++)
				std::cout<< par[i] << Form("  par%d \n", i);
			IntSkewNorm = 1.0;
		}
		double xi = par[6];
		double sigma = par[4];
		double skew  = par[5];
		boost::math::skew_normal dist(xi, sigma, skew);
	 	IntSkewNorm = boost::math::cdf(dist, xmax) - boost::math::cdf(dist, xmin);
*/

		double alpha = abs(par[2]+1.0e-6);
		double n = par[3];
		double mean = par[1];
		double sigma = par[0];

		IntCB = TMath::Abs(-ROOT::Math::crystalball_integral(2.0*mean-min, alpha, n, sigma, mean) + ROOT::Math::crystalball_integral(2.0*mean-max, alpha, n, sigma, mean));
	}
}

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
