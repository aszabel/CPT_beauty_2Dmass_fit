#include "M_B_2missPT_fit.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

namespace cpt_b0_analysis{

   mb_2misspt_fit::mb_2misspt_fit(double _minn, double _maxx, int _ncontr):minn(_minn), maxx(_maxx), ncontr(_ncontr){ 
	   for (int i=0; i<ncontr; i++){
		   int_cos.push_back(1.0);
		   int_gaus1.push_back(1.0);
		   int_gaus2.push_back(1.0);
	   }
	}
double mb_2misspt_fit::func_full(const double *xx, const double *par, int icontr){
                   auto raised_cosine = [this, icontr](const double *x, const double *par)->double{
                        double mass = x[0];
                        double s = abs(par[1]);
                        double mean = par[0];
                        double cos = 0.0;
                        if (mass>=mean-s && mass<=mean+s){
                                cos = 1.0/(2.0*s)*(1.0+TMath::Cos((mass-mean)/s*TMath::Pi()));
                        }
                        double minnB = mean-s;
                        if (minn>mean-s) minnB = minn;
                        double maxxB = mean+s;
                        if (maxx<mean+s) maxxB = maxx;
                        //double int_cos = 1.0/(2.0)*(1.0+(maxxB-mean)/s+TMath::Sin((maxxB-mean)/s*TMath::Pi())/TMath::Pi())-1.0/(2.0)*(1.0+(minnB-mean)/s+TMath::Sin((minnB-mean)/s*TMath::Pi())/TMath::Pi());
                	if (int_cos[icontr]!=0.0) cos /=int_cos[icontr];

                        return  cos;
                };
		auto gaus = [this, icontr](const double *x, const double *par)->double{
                        const double m_rec = x[0];
                        double norm1 = 1.0 -abs(par[2]);
                        double mean1 = par[3];
                        double sigma1 = abs(par[4]);
			//cout << sigma1 << "<==sigma1  mean1==> " << mean1 << endl;
                        double gaus1 = ROOT::Math::gaussian_pdf(m_rec, sigma1, mean1);
                        //double int_gaus1 =  ROOT::Math::normal_cdf(maxx, sigma1, mean1)-ROOT::Math::normal_cdf(minn, sigma1, mean1);
			gaus1*=abs(norm1);
                        if(int_gaus1[icontr]!=0) gaus1/=int_gaus1[icontr];
                        double norm2 = abs(par[2]);
                        double mean2 = par[5];
                        double sigma2 = abs(par[6]);
			//cout << sigma2 << "<==sigma2  mean2==> " << mean2 << endl;
                        double gaus2 = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean2);
                        //double int_gaus2 =  ROOT::Math::normal_cdf(maxx, sigma2, mean2)-ROOT::Math::normal_cdf(minn, sigma2, mean2);
			gaus2*=abs(norm2);
			if(int_gaus2[icontr]!=0.0) gaus2/=int_gaus2[icontr];
                        return (gaus1+gaus2);
                };

	return (raised_cosine(xx, par)+gaus(xx, par))/2.0;
}
}
