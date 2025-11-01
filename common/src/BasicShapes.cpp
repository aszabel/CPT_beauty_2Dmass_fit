#include "BasicShapes.h"

#include <Math/DistFunc.h>	// normal_cdf

#include <boost/math/special_functions/owens_t.hpp>

#include "Math/GaussIntegrator.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/WrappedTF1.h"
#include "TF1.h"
#include "TMath.h"

namespace cpt_b0_analysis {

double RaisedCosinePDF::EvalPDF(const double *xx, const double *par, const int component) {
	double mass = xx[0];
	double s = abs(par[1]);
	double mean = par[0];
	double cos = 0.0;
	if (mass >= mean - s && mass <= mean + s) {
		cos = 1.0 / (2.0 * s) * (1.0 + TMath::Cos((mass - mean) / s * TMath::Pi()));
	}
	if (IntCos != 0.0) cos /= IntCos;
	return cos;
}

void RaisedCosinePDF::CalcIntegral(const double *par, double min, double max) {
	double s = abs(par[1]);
	double mean = par[0];
	double minnB = mean - s;
	if (min > mean - s) minnB = min;
	double maxxB = mean + s;
	if (max < mean + s) maxxB = max;

	IntCos =
		1.0 / (2.0) *
			(1.0 + (maxxB - mean) / s +
			 TMath::Sin((maxxB - mean) / s * TMath::Pi()) / TMath::Pi()) -
		1.0 / (2.0) *
			(1.0 + (minnB - mean) / s + TMath::Sin((minnB - mean) / s * TMath::Pi()) / TMath::Pi());
}

double GaussPDF::EvalPDF(const double *xx, const double *par, const int component) {
	const double m_rec = xx[0];
	double mean = par[0];
	double sigma = abs(par[1]);
	double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
	if (IntGaus != 0.0) {
		gaus /= IntGaus;
	}

	return gaus;
}

void GaussPDF::CalcIntegral(const double *par, double min, double max) {
	double mean = par[0];
	double sigma = abs(par[1]);

	IntGaus = ROOT::Math::normal_cdf(max, sigma, mean) - ROOT::Math::normal_cdf(min, sigma, mean);
}

/*
	double DoubleGaussPDF::EvalPDF(const double *xx, const double *par, const
   int component)
	{
		auto gaus = [this](const double *x, const double *par, const int idx) ->
   double
		{
			const double m_rec = x[0];
			double norm = 0.0;
			if (idx == 0)
				norm = 1.0 - abs(par[0]);
			else
				norm = abs(par[0]);
			if (norm<=1.0e-9)
				norm = 0.0;

			double mean = abs(par[1 + idx * 2]);
			double sigma = abs(par[2 + idx * 2]);
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

		if (component == 0)
			return gaus(xx, par, 0);
		else if (component == 2)
			return gaus(xx, par, 1);
		else
			return gaus(xx, par, 0) + gaus(xx, par, 1);
	}

	void DoubleGaussPDF::CalcIntegral(const double *par, double min, double max)
	{
		double mean1 = abs(par[3]);
		double sigma1 = abs(par[4]);

		IntGaus1 = ROOT::Math::normal_cdf(max, sigma1, mean1) -
   ROOT::Math::normal_cdf(min, sigma1, mean1);

		if (abs(par[2])-1.0 < 1.0e-6){
			IntGaus1 = 0.0;
		}
		double mean2 = abs(par[5]);
		double sigma2 = abs(par[6]);

		IntGaus2 = ROOT::Math::normal_cdf(max, sigma2, mean2) -
   ROOT::Math::normal_cdf(min, sigma2, mean2); if (abs(par[2])<1.0e-6) IntGaus2
   = 0.0;
	}
*/

double DoubleGaussPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double f = abs(par[0]);
	double par1[2];
	double par2[2];
	par1[0] = par[1];
	par1[1] = par[2];
	par2[0] = par[3];
	par2[1] = par[4];

	if (component == 0)
		return (1 - f) * Gauss1.EvalPDF(xx, par1);
	else if (component == 1)
		return f * Gauss2.EvalPDF(xx, par2);
	else
		return (1 - f) * Gauss1.EvalPDF(xx, par1) + f * Gauss2.EvalPDF(xx, par2);
}

void DoubleGaussPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[2];
	double par2[2];
	par1[0] = par[1];
	par1[1] = par[2];
	par2[0] = par[3];
	par2[1] = par[4];

	Gauss1.CalcIntegral(par1, min, max);
	Gauss2.CalcIntegral(par2, min, max);
}

double SkewNormalPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double xm = xx[0] - par[0];	 // mean
	double sigma = abs(par[1]);
	double skew = par[2];
	double sn = 0.0;

	if (skew == 0.0)
		sn = ROOT::Math::gaussian_pdf(xm, sigma, 0.0);
	else
		sn = 2. * ROOT::Math::gaussian_pdf(xm, sigma, 0.0) *
			 ROOT::Math::normal_cdf(skew * xm / sigma, 1.0, 0.0);

	if (IntSkewNorm != 0.0) {
		sn /= IntSkewNorm;
	}

	return sn;
}

void SkewNormalPDF::CalcIntegral(const double *par, double min, double max) {
	double xmin = min - par[0];
	double xmax = max - par[0];
	double sigma_sk = abs(par[1]);
	double skew = par[2];

	double zmin = xmin / sigma_sk;
	double zmax = xmax / sigma_sk;
	double Phi = ROOT::Math::normal_cdf(zmax) - ROOT::Math::normal_cdf(zmin);

	if (skew == 0) {
		IntSkewNorm = Phi;
	} else {
		double T = boost::math::owens_t(zmax, skew) - boost::math::owens_t(zmin, skew);
		IntSkewNorm = Phi - 2.0 * T;
	}
}

double CrystalBallPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double m_rec = xx[0];
	double mean = par[0];
	double sigma = abs(par[1]);
	double alpha = abs(par[2]) + 1.0e-10;
	double n = par[3];

	double CB = ROOT::Math::crystalball_function(2.0 * mean - m_rec, alpha, n, sigma, mean);
	if (IntCB != 0.0) CB /= IntCB;
	return CB;
}

void CrystalBallPDF::CalcIntegral(const double *par, double min, double max) {
	double mean = par[0];
	double sigma = abs(par[1]);
	double alpha = abs(par[2]) + 1.0e-10;
	double n = par[3];

	IntCB = TMath::Abs(-ROOT::Math::crystalball_integral(2.0 * mean - min, alpha, n, sigma, mean) +
					   ROOT::Math::crystalball_integral(2.0 * mean - max, alpha, n, sigma, mean));
}

double JohnsonPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double xi = par[0];
	double lambda = abs(par[1]);
	double gamma = par[2];
	double delta = abs(par[3]);

	double z = (xx[0] - xi) / lambda;
	double asinh_z = TMath::ASinH(z);
	double arg = gamma + delta * asinh_z;

	double norm = delta / (lambda * TMath::Sqrt(TMath::TwoPi()));
	double denom = TMath::Sqrt(1.0 + z * z);
	double expo = TMath::Exp(-0.5 * arg * arg);

	double JSU = norm / denom * expo;
	if (IntJSU != 0.0) JSU /= IntJSU;
	return JSU;
}

void JohnsonPDF::CalcIntegral(const double *par, double min, double max) {
	double xi = par[0];
	double lambda = abs(par[1]);
	double gamma = par[2];
	double delta = abs(par[3]);
	double argmax = gamma + delta * TMath::ASinH((max - xi) / lambda);
	double argmin = gamma + delta * TMath::ASinH((min - xi) / lambda);
	IntJSU =
		abs(ROOT::Math::normal_cdf(argmax, 1.0, 0.0) - ROOT::Math::normal_cdf(argmin, 1.0, 0.0));
}

double DoubleSidedCrystalBallPDF::EvalPDF(const double *xx, const double *par,
										  const int component) {
	double mean = par[0];
	double sigma = abs(par[1]);
	double alpha = abs(par[2]);
	double n = par[3];
	double alpha_h = abs(par[4]);
	double n2 = par[5];
	double m_rec = xx[0];
	double result;

	if (m_rec < mean) {
		result = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
	} else {
		result = ROOT::Math::crystalball_function(2. * mean - m_rec, alpha_h, n2, sigma, mean);
	}

	if (IntDCB != 0) result /= IntDCB;
	return result;
}

void DoubleSidedCrystalBallPDF::CalcIntegral(const double *par, double min, double max) {
	double alpha = par[2];
	double n = par[3];
	double mean = par[0];
	double sigma = abs(par[1]);
	double alpha_h = abs(par[4]);
	double n2 = par[5];

	IntDCB =
		TMath::Abs(-ROOT::Math::crystalball_integral(min, alpha, n, sigma, mean) +
				   ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) +
		TMath::Abs(-ROOT::Math::crystalball_integral(2. * mean - max, alpha_h, n2, sigma, mean) +
				   ROOT::Math::crystalball_integral(mean, alpha_h, n2, sigma, mean));
}
}  // namespace cpt_b0_analysis

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
