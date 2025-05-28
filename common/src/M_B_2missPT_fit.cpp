#include "M_B_2missPT_fit.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis
{

	RaisedCosinePlusGaussPDF::RaisedCosinePlusGaussPDF()
	{
		IntCos = 1.0;
		IntGaus1 = 1.0;
		IntGaus2 = 1.0;
	}

	double RaisedCosinePlusGaussPDF::EvalPDF(const double *xx, const double *par)
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
		auto gaus = [this](const double *x, const double *par) -> double
		{
			const double m_rec = x[0];
			double norm1 = 1.0 - abs(par[2]);
			if (norm1<=1.0e-6) 
				norm1 = 0.0;
			double mean1 = abs(par[3]);
			double sigma1 = abs(par[4]);
			double gaus1 = ROOT::Math::gaussian_pdf(m_rec, sigma1, mean1);
			gaus1 *= norm1;
			if (abs(IntGaus1) >= 1.0e-6)
				gaus1 /= IntGaus1;
			double norm2 = abs(par[2]);
			if (norm2<=1.0e-9) 
				norm2 = 0.0;
			double mean2 = abs(par[5]);
			double sigma2 = abs(par[6]);
			double gaus2 = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean2);
			gaus2 *= norm2;
			if (IntGaus2 != 0.0)
				gaus2 /= IntGaus2;
			return (gaus1 + gaus2);
		};
		return (raised_cosine(xx, par) + gaus(xx, par)) / 2.0;
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

	SkewNormalPlusGausPDF::SkewNormalPlusGausPDF()
	{
		IntSkewNorm = 1.0;
		IntCB = 1.0;
	       	skew_normal = [this](const double *x, const double *par) -> double
                {
			double xm = x[0]-par[6];
		        double sigma = par[4];
        		double skew  = par[5];

			double skew_norm = 2.*ROOT::Math::gaussian_pdf(xm, sigma, 0.0)*ROOT::Math::normal_cdf(skew*xm/sigma,1.0, 0.0);
                        return skew_norm;
        	};


	}

        double SkewNormalPlusGausPDF::EvalPDF(const double *xx, const double *par)
        {
                auto gaussian = [this](const double *x, const double *par) -> double
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
		return abs(1.0-par[7])*skew_normal(xx, par)/IntSkewNorm+abs(par[7])*gaussian(xx, par);
	}
       void SkewNormalPlusGausPDF::CalcIntegral(const double *par, double min, double max)
       {
	       
	       TF1 skewfunc("skewfunc", skew_normal, min, max, 8);
               skewfunc.SetParameters(par);
               ROOT::Math::WrappedTF1 wf1(skewfunc);
               ROOT::Math::GaussIntegrator ig;
               ig.SetFunction(wf1);
               ig.SetRelTolerance(0.01);
               IntSkewNorm = ig.Integral(min, max);
	

		
       
        
                double alpha = abs(par[2]+1.0e-6);
                double n = par[3];
                double mean = par[1];
                double sigma = par[0];

                IntCB = TMath::Abs(-ROOT::Math::crystalball_integral(2.0*mean-min, alpha, n, sigma, mean) + ROOT::Math::crystalball_integral(2.0*mean-max, alpha, n, sigma, mean));
      }

}
