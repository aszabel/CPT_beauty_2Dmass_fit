#include "M_B_2missPT_fit.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis{

RaisedCosinePlusGaussPDF::RaisedCosinePlusGaussPDF(){ 
		   IntCos = 1.0;
		   IntGaus1 = 1.0;
		   IntGaus2 = 1.0;
}
double RaisedCosinePlusGaussPDF::EvalPDF(const double *xx, const double *par){
                   auto raised_cosine = [this](const double *x, const double *par)->double{
                        double mass = x[0];
                        double s = abs(par[1]);
                        double mean = par[0];
                        double cos = 0.0;
                        if (mass>=mean-s && mass<=mean+s){
                                cos = 1.0/(2.0*s)*(1.0+TMath::Cos((mass-mean)/s*TMath::Pi()));
                        }
        		if (IntCos !=0.0) cos /=IntCos;
                        return  cos;
                };
		auto gaus = [this](const double *x, const double *par)->double{
                        const double m_rec = x[0];
                        double norm1 = 1.0 -abs(par[2]);
                        double mean1 = par[3];
                        double sigma1 = abs(par[4]);
			double gaus1 = ROOT::Math::gaussian_pdf(m_rec, sigma1, mean1);
			gaus1*=abs(norm1);
                        if(IntGaus1!=0) gaus1/=IntGaus1;
                        double norm2 = abs(par[2]);
                        double mean2 = par[5];
                        double sigma2 = abs(par[6]);
			//cout << sigma2 << "<==sigma2  mean2==> " << mean2 << endl;
                        double gaus2 = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean2);
			gaus2*=abs(norm2);
			if(IntGaus2!=0.0) gaus2/=IntGaus2;
                        return (gaus1+gaus2);
                };

	return (raised_cosine(xx, par)+gaus(xx, par))/2.0;
}
void RaisedCosinePlusGaussPDF::CalcIntegral(const double *par, double min, double max){

        double s = abs(par[1]);
        double mean = par[0];
        double minnB = mean-s;
        if (min>mean-s) minnB = min;
        double maxxB = mean+s;
        if (max<mean+s) maxxB = max;

	IntCos = 1.0/(2.0)*(1.0+(maxxB-mean)/s+TMath::Sin((maxxB-mean)/s*TMath::Pi())/TMath::Pi())-1.0/(2.0)*(1.0+(minnB-mean)/s+TMath::Sin((minnB-mean)/s*TMath::Pi())/TMath::Pi());

        double mean1 = par[3];
        double sigma1 = abs(par[4]);

        IntGaus1 = ROOT::Math::normal_cdf(max, sigma1, mean1)-ROOT::Math::normal_cdf(min, sigma1, mean1);

        double mean2 = par[5];
        double sigma2 = abs(par[6]);

        IntGaus2 = ROOT::Math::normal_cdf(max, sigma2, mean2)-ROOT::Math::normal_cdf(min, sigma2, mean2);

}

}
