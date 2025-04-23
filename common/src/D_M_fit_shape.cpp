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

	DoubleSidedCrystalballPlusGaussPDF::DoubleSidedCrystalballPlusGaussPDF()
	{
		IntGaus = 1.0;
		IntDCB = 1.0;
	}
	double DoubleSidedCrystalballPlusGaussPDF::EvalPDF(const double *xx, const double *_par)
	{
		auto gausPDF = [this](const double *x, const double *par) -> double
		{
			double m_rec = x[0];
			double sigma = par[0];
			double mean = par[1];
			double f12 = abs(par[3]);
			double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
			//cout<<"Gaus: "<<gaus<<endl;
			//cout<<"IntGaus: "<<IntGaus<<endl;
			if (IntGaus != 0)
				gaus /= IntGaus;

			return f12 * gaus;
		};

		auto DoubleSidedCrystalballFunction = [this](const double *x, const double *par) -> double
		{
			double alpha = par[4];
			double n = par[5];
			double n2 = par[7];
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
				result = ROOT::Math::crystalball_function(2. * mean - m_rec, alpha_h, n2, sigma, mean);
			}

			//cout<<"DCB: "<<result<<endl;
			//cout<<"IntDCB: "<<IntDCB<<endl;
			if (IntDCB != 0)
				result /= IntDCB;
			return (1.0 - f12) * result;
		};

		return gausPDF(xx, _par) + DoubleSidedCrystalballFunction(xx, _par);
	}
	void DoubleSidedCrystalballPlusGaussPDF::CalcIntegral(const double *par, double min, double max)
	{
		double sigma = par[0];
		double mean = par[1];

		IntGaus = ROOT::Math::normal_cdf(max, sigma, mean) - ROOT::Math::normal_cdf(min, sigma, mean);

		double alpha = par[4];
		double n = par[5];
		mean = par[1];
		sigma = par[2];
		double alpha_h = abs(par[6]);
		double n2 = par[7];

		IntDCB = TMath::Abs(-ROOT::Math::crystalball_integral(min, alpha, n, sigma, mean) + ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2. * mean - max, alpha_h, n2, sigma, mean) + ROOT::Math::crystalball_integral(mean, alpha_h, n2, sigma, mean));
	}

        JohnsonPlusGaussPDF::JohnsonPlusGaussPDF()
        {
                IntJSU = 1.0;
                IntGaus = 1.0;
        }
        double JohnsonPlusGaussPDF::EvalPDF(const double *xx, const double *_par)
        {
		auto gausPDF = [this](const double *x, const double *par) -> double
                {
                        double m_rec = x[0];
                        double sigma = par[0];
                        double mean = par[1];
                        double f12 = abs(par[3]);
                        double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
                        //cout<<"Gaus: "<<gaus<<endl;
                        //cout<<"IntGaus: "<<IntGaus<<endl;
                        if (IntGaus != 0)
                                gaus /= IntGaus;

                        return f12 * gaus;
		};

		auto JohnsonPDF = [this](const double *x, const double *par) -> double
		{

			double gamma = par[4]; 
 			double delta = par[5]; 
  			double xi    = par[1]; 
  			double lambda= abs(par[2]);
			double f12 = abs(par[3]);

  			double z = (x[0] - xi) / lambda;
  			double asinh_z = TMath::ASinH(z);
 			double arg = gamma + delta * asinh_z;

  			double norm = delta / (lambda * TMath::Sqrt(TMath::TwoPi()));
  			double denom = TMath::Sqrt(1.0 + z * z);
  			double expo = TMath::Exp(-0.5 * arg * arg);

  			//return (1-f12) * norm * (1.0 / denom) * expo;
			
  			double JSU = norm*  (1-f12) / denom * expo;
			if (IntJSU!=0.0)
				JSU/= IntJSU;
			return JSU;
		};
		return gausPDF(xx, _par)+JohnsonPDF(xx, _par);

	}
	void JohnsonPlusGaussPDF::CalcIntegral(const double *par, double min, double max)
	{
		double sigma = par[0];
		double mean = par[1];

		IntGaus = ROOT::Math::normal_cdf(max, sigma, mean) - ROOT::Math::normal_cdf(min, sigma, mean);

                double gamma = par[4];
                double delta = par[5];
                double xi    = par[1];
                double lambda= abs(par[2]);
		double argmax = gamma + delta * TMath::ASinH((max - xi) / lambda);
		double argmin = gamma + delta * TMath::ASinH((min - xi) / lambda);
                //IntJSU = 0.5 * (1.0 + TMath::Erf(zmax / TMath::Sqrt2())) - 0.5 * (1.0 + TMath::Erf(zmin/ TMath::Sqrt2()));
		IntJSU = abs(ROOT::Math::normal_cdf(argmax, 1.0, 0.0) - ROOT::Math::normal_cdf(argmin, 1.0, 0.0));
	}

}
