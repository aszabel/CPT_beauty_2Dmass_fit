#include "D_M_fit_shape.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis{

md_fit::md_fit(double _minx, double _maxx):minx(_minx), maxx(_maxx){}

double md_fit::func_full(const double *xx, const double *_par){
	auto gausPDF = [this](const double *x, const double *par)->double{
  		double m_rec = x[0];
  		double sigma = par[0];
  		double mean = par[1];
  		double f12 = abs(par[3]);
  		double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  		double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean);
  		gaus/=int_gaus;

  		return f12*gaus;
	};

	auto DoubleSidedCrystalballFunction = [this](const double *x, const double *par)->double{
  		double alpha = par[4];
  		double n     = par[5];
  		double mean  = par[1];
  		double sigma = par[2];
  		double f12 = abs(par[3]);
  		double alpha_h = abs(par[6]);
  		double m_rec=x[0];
  		double result;

   		if (m_rec<mean){
       			result = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
   		}
   		else{
       			result = ROOT::Math::crystalball_function(2.*mean-m_rec, alpha_h, n, sigma, mean);
   		}
      			double intDCB = TMath::Abs(-ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2.*mean-maxx, alpha_h, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha_h, n, sigma, mean));

   		return (1.0-f12)*result/intDCB;
		};
	return gausPDF(xx, _par)+DoubleSidedCrystalballFunction(xx, _par);
}
}
