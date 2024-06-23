#include "ChebyshevPDF.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis{

ChebyshevPDF::ChebyshevPDF(){
		IntCheb = 1.0;
}
double ChebyshevPDF::EvalPDF(const double *xx, const double *_par){
	auto DebPDF = [this](const double *x, const double *par)->double{
		double m_rec = x[0];
 	        double a1 = par[0];
  		double a2 = par[1];
  		double cheb = 1.0 + a1*m_rec + a2*(2.0*m_rec*m_rec-1.0);
  		cheb/=IntCheb;
  		return cheb;
	};

	return DebPDF(xx, _par);
}
void ChebyshevPDF::CalcIntegral(const double *par, double min, double max){
        double a1 = par[0];
        double a2 = par[1];
 	
	IntCheb = (1.0-a2)*max + 0.5*a1*max*max+2./3.*a2*max*max*max- (1.0-a2)*min-0.5*a1*min*min-2./3.*a2*min*min*min;
}

}
