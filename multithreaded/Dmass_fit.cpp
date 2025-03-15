/**
  * Adam Szabelski 22.01.2020
  * 
  * @Dmass_fit.cpp - implementation of fit invariant mass of D+(-)-> K pi pi. The fraction of combinatorial background is a common parameter
  *                  with the lifetime description
**/
#include "Dmass_fit.h"

#include "TMath.h"
#include "Math/Math.h"
#include "Math/Functor.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TStyle.h"
#include "TFile.h"

#include <iostream>

namespace cpt_b0_analysis{


double gaus_pdf (double *x, double *par){
  int nbins = 20;
  double bin_width = 12./int(sqrt(nbins));
  double sigma = par[0];
  double mean = par[1];
  return double(nbins)*bin_width*ROOT::Math::gaussian_pdf(x[0], sigma, mean);
}

Dmass::Dmass(std::vector<std::pair<double, double>>& vect_Dmass_, TH2D **hcorr_,  double minx_, double maxx_, double miny_, double maxy_, bool sign_) : vect_Dmass(vect_Dmass_), hcorr(hcorr_), minx(minx_), maxx(maxx_), miny(miny_), maxy(maxy_), sign(sign_){
   nbins = 40;
   bin_width = (maxx-minx)/double(nbins);
   bin_widthy = (maxy-miny)/double(nbins);
   nevents = vect_Dmass.size();
   hdata.SetBins(nbins, minx, maxx, nbins, miny, maxy);
   hdata.SetName("hdata");
   for (const auto& ivect: vect_Dmass){
	   hdata.Fill(get<0>(ivect), get<1>(ivect));
   }
   cout << hdata.Integral() << " <++++++++++++++++++++ tuuuuu\n";

   mb_2misspt_graphs = new mb_2misspt_fit( miny, maxy);
}

Dmass::~Dmass(){}

void Dmass::SetExternal(double xx_[ncontr][nprm], double dxx_[ncontr][nprm], double yy_[nsigBu+nBs+nB0D+nBuD], double dyy_[nsigBu+nBs+nB0D+nBuD]){
	for (int i=0; i<ncontr; i++){
	  for (int j=0; j<nprm; j++){
		xx[i][j] = xx_[i][j]; 
		dxx[i][j] = dxx_[i][j];
		//cout << xx[i][j] << "  tam tam " << dxx[i][j] << endl;
	  }
	}
	for (int i=0; i < nsigBu+nBs+nB0D+nBuD; i++){
        	yy[i] = yy_[i];	
        	dyy[i] = dyy_[i];	
		//cout << yy[i] << "  tamyy tamyy " << dyy[i]<< endl;
	}
}	
double Dmass::operator()(const double *par){
  //cout << " OP \n";
  double chi2 = calcChi2(par); 
  //cout << chi2<< " chi2_result " << endl;
  return chi2;

}
double Dmass::gausPDF(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = abs(par[0]);
  double mean = par[1];
  double f12 = abs(par[2]);
  //double nevents = par[24];
  double frac_sig = 1.0-abs(par[82]);
  double frac_Bu = abs(par[82])*(1.0-abs(par[83]));
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean);
  gaus/=int_gaus;
  TH2D hist(*hcorr[0]);
  //hist.Scale(1.0-per_Bu_sig);
  //hist.Add(hcorr[1], per_Bu_sig);
  //hist.Scale(1./hist.Integral());
  //normalize_histogram(hist);
  //cout << mb_2misspt_graphs->func_full(m_corr, &par[28], 0) << endl; 
  //cout << mb_2misspt_graphs->func_full(m_corr, &par[28], 0) << " <========> " << mb_2misspt_graphs->func_full(m_corr, &par[28], 1)<<endl;
  //cout << abs(1.0-fcomb) <<  "  " <<  abs(1.0-fD-fDCB-fCB)  <<"  ddddddddddddd " << f12 << endl ;
  double cont  = frac_sig*mb_2misspt_graphs->func_full(&x[1], &par[28], 0) + frac_Bu*mb_2misspt_graphs->func_full(&x[1], &par[28], 1);
  //cout << cont << " contgaus\n";
  return bin_width*double(nevents)*f12*gaus*cont;
}

double Dmass::crystalballPDF(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double mean = par[1];
  double sigma2 = abs(par[3]);
  double f12 = abs(par[2]);
  double alpha = abs(par[4]);
  double n = par[5];
  double frac_sig = 1.0-abs(par[82]);
  double frac_Bu = abs(par[82])*(1.0-abs(par[83]));
  //double nevents = par[18];
  //double CB = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean);
  //double int_CB =  ROOT::Math::normal_cdf(maxx, sigma2, mean)-ROOT::Math::normal_cdf(minx, sigma2, mean);
  double CB = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma2, mean);
  double int_CB = TMath::Abs(ROOT::Math::crystalball_integral(minx, alpha, n, sigma2, mean)-ROOT::Math::crystalball_integral(maxx, alpha, n, sigma2, mean));
  CB/=int_CB;

  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean);
  double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma2, mean)-ROOT::Math::normal_cdf(minx, sigma2, mean);
  gaus/=int_gaus;
  
  TH2D hist(*hcorr[0]);
   //double per_Bu_sig = frac_Bu/(frac_Bu+frac_sig);
  //hist.Scale(1.0-per_Bu_sig);
  //hist.Add(hcorr[1], per_Bu_sig);
  //hist.Scale(1./hist.Integral());
  //normalize_histogram(hist);

  //double cont  = (1.-per_Bu_sig)*mb_2misspt_graphs->func_full(&x[1], &par[28], 0) + per_Bu_sig*mb_2misspt_graphs->func_full(&x[1], &par[28], 1);
  double cont  = frac_sig*mb_2misspt_graphs->func_full(&x[1], &par[28], 0) + frac_Bu*mb_2misspt_graphs->func_full(&x[1], &par[28], 1);
  //cout << cont << " CB \n";
  //return bin_width*double(nevents)*abs(1.0-fcomb)*abs(1.0-fD-fCB-fDCB)*abs(1.0-f12)*gaus*cont;
  return bin_width*double(nevents)*(1.0-f12)*CB*cont;
}
double Dmass::DoubleSidedCrystalballFunction_BuD(double *x, const double *par, int binx, int bin_y)
{
  double alpha = abs(par[21]);
  double n     = par[22];
  double mean  = par[1];
  double sigma = abs(par[19]);
  double f12 = abs(par[20]);
  double alpha_h = abs(par[21]);
  double m_rec=x[0];
  double m_corr = x[1];
  double result;
  double frac_DuB = abs(par[82])*abs(par[83])*abs(par[84])*abs(par[85])*abs(par[86]);
   if (m_rec<mean){
       result = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
   }
   else{
       result = ROOT::Math::crystalball_function(2.*mean-m_rec, alpha_h, n, sigma, mean);
   }
      double intDCB = TMath::Abs(-ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2.*mean-maxx, alpha_h, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha_h, n, sigma, mean));
  
  //double cont  = hcorr[3]->GetBinContent(bin_x, bin_y);
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 5);
 
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean);
 
  gaus/=int_gaus;


  double CB = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
  double int_CB = TMath::Abs(ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean)-ROOT::Math::crystalball_integral(maxx, alpha, n, sigma, mean));
  CB/=int_CB;
   return bin_width*double(nevents)*frac_DuB*(1.-f12)*result/intDCB*cont;
   //return bin_width*double(nevents)*(1.0-fcomb)*fDCB*result/intDCB*cont;
}


double Dmass::gausPDF_BuD(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = par[17];
  double mean = par[1];
  double f12 = abs(par[20]);
  double frac_DuB = abs(par[82])*abs(par[83])*abs(par[84])*abs(par[85])*abs(par[86]);
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean);
  gaus/=int_gaus;
  //double cont  = hcorr[5]->GetBinContent(bin_x, bin_y);
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 5);

  //cout << cont << " DspDm \n";
  return bin_width*double(nevents)*frac_DuB*f12*gaus*cont;
}


double Dmass::gausPDF_B0D(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = abs(par[10]);
  double mean = par[1];
  double f12 = abs(par[13]);
  double frac_B0D = abs(par[82])*abs(par[83])*abs(par[84])*(1.0-abs(par[85]));
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean);
  gaus/=int_gaus;
  //double cont  = hcorr[5]->GetBinContent(bin_x, bin_y);
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 3);

  //cout << cont << " DspDm \n";
  return bin_width*double(nevents)*frac_B0D*f12*gaus*cont;
}
double Dmass::DoubleSidedCrystalballFunction_B0D(double *x, const double *par, int binx, int bin_y)
{
  double alpha = abs(par[14]);
  double n     = par[15];
  double mean  = par[1];
  double sigma = abs(par[12]);
  double f12 = abs(par[13]);
  double alpha_h = abs(par[14]);
  double m_rec=x[0];
  double m_corr = x[1];
  double result;
  double frac_B0D = abs(par[82])*abs(par[83])*abs(par[84])*(1.0-abs(par[85]));
   if (m_rec<mean){
       result = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
   }
   else{
       result = ROOT::Math::crystalball_function(2.*mean-m_rec, alpha_h, n, sigma, mean);
   }
      double intDCB = TMath::Abs(-ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2.*mean-maxx, alpha_h, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha_h, n, sigma, mean));
  //cout << intDCB << "  int DCB " << result << endl; 
  //double cont  = hcorr[3]->GetBinContent(bin_x, bin_y);
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 3);
  
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(maxx, sigma, mean)-ROOT::Math::normal_cdf(minx, sigma, mean);
 
  gaus/=int_gaus;


  double CB = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
  double int_CB = TMath::Abs(ROOT::Math::crystalball_integral(minx, alpha, n, sigma, mean)-ROOT::Math::crystalball_integral(maxx, alpha, n, sigma, mean));
  CB/=int_CB;
   return bin_width*double(nevents)*frac_B0D*(1.0-f12)*result/intDCB*cont;
   //return bin_width*double(nevents)*(1.0-fcomb)*fDCB*result/intDCB*cont;
}

double Dmass::crystalballPDF2(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double mean = par[7];
  double sigma2 = par[6];
  double alpha = par[8];
  double n = par[9];
  double frac_Bs = abs(par[82])*abs(par[83])*(1.0-abs(par[84]));
  //double CB = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean);

  //double int_CB =  ROOT::Math::normal_cdf(maxx, sigma2, mean)-ROOT::Math::normal_cdf(minx, sigma2, mean);
  double CB = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma2, mean);
  double int_CB = TMath::Abs(ROOT::Math::crystalball_integral(minx, alpha, n, sigma2, mean)-ROOT::Math::crystalball_integral(maxx, alpha, n, sigma2, mean));
  CB/=int_CB;
  //return bin_width*double(nevents)*(1.0-fcomb)*(1.0-f12)*abs(1.0-fCB)*abs(1.0-fDCB)*CB;
  
  //double cont  = hcorr[2]->GetBinContent(bin_x, bin_y);
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 2);
  return bin_width*double(nevents)*frac_Bs*CB*cont;
}
double Dmass::exponentPDF(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double exp_slope = par[24];
  double frac_exp = abs(par[82])*abs(par[83])*abs(par[84])*abs(par[85])*(1.0-abs(par[86]));  //double expo = TMath::Exp(-exp_slope*m_rec);
  //double int_expo = 1.0/exp_slope *(TMath::Exp(-exp_slope*minx)-TMath::Exp(-exp_slope*maxx));
  double expo = ROOT::Math::exponential_pdf(m_rec, exp_slope);
  double int_expo =  ROOT::Math::exponential_cdf(maxx, exp_slope)-ROOT::Math::exponential_cdf(minx, exp_slope);
  expo /= int_expo;
 
  //double cont  = hcorr[4]->GetBinContent(bin_x, bin_y);
  //double cont  = tabcorr[4][bin_x-1][bin_y-1];
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 4);
  return bin_width*double(nevents)*frac_exp*expo*cont;
}
double Dmass::chebyshevPDF(double *x, const double *par, int bin_x, int bin_y){
  double m_rec = x[0];
  double m_corr = x[1];
  double a1 = par[25];
  double a2 = par[26];
  double frac_cheb = abs(par[82])*abs(par[83])*abs(par[84])*abs(par[85])*(1.0-abs(par[86]));
  //double nevents = par[24];
  double cheb = 1.0 + a1*m_rec + a2*(2.0*m_rec*m_rec-1.0);
  double int_cheb = (1.0-a2)*maxx + 0.5*a1*maxx*maxx+2./3.*a2*maxx*maxx*maxx- (1.0-a2)*minx-0.5*a1*minx*minx-2./3.*a2*minx*minx*minx;
  cheb/=int_cheb;
  
  //double cont  = hcorr[4]->GetBinContent(bin_x, bin_y);
  //double cont  = tabcorr[4][bin_x-1][bin_y-1];
  double cont  = mb_2misspt_graphs->func_full(&x[1], &par[28], 4);
  //cout << cont << " cheb cont \n";
  return bin_width*double(nevents)*frac_cheb*cheb*cont;
}
void Dmass::GetCorrBin(int &bin_x, int &bin_y, double m_rec, double m_corr){
  double Nbins_x = (double)hcorr[0]->GetNbinsX();  
  double Nbins_y = (double)hcorr[0]->GetNbinsY();

  
  double bin_width_corr_x = (maxx-minx)/Nbins_x;
  bin_x = int((m_rec-minx)/bin_width_corr_x+1); 

  double bin_width_corr_y = (maxy-miny)/Nbins_y;
  bin_y = int((m_corr-miny)/bin_width_corr_y+1);
} 

double Dmass::Dmass_func(double *x, const double *par){

  double m_rec = x[0];
  double m_corr = x[1];

  int bin_x, bin_y;
  GetCorrBin(bin_x, bin_y, m_rec, m_corr); 

  double expo = exponentPDF(x, par, bin_x, bin_y);
  //cout << expo << " expo\n";
  double gaus1 = gausPDF(x, par, bin_x, bin_y);
  //cout << gaus1 << " gaus1\n";
  double gausB0D = gausPDF_B0D(x, par, bin_x, bin_y);
  //cout << gausB0D << " gausD0B\n";
  double gausBuD = gausPDF_BuD(x, par, bin_x, bin_y);
  //cout << gausBuD << " gausDuB\n";
  double gaus2 =  crystalballPDF(x, par, bin_x, bin_y);
  //cout << gaus2 << " gaus2\n";
  double CBs = crystalballPDF2(x, par, bin_x, bin_y);
  //cout << CBs << " CBs\n";
  double cheb = chebyshevPDF(x, par, bin_x, bin_y);
  //cout << cheb << " cheb\n";
  double DCB_B0D = DoubleSidedCrystalballFunction_B0D(x, par, bin_x, bin_y);
  //cout << DCB_B0D << " DCB_D0B\n";
  double DCB_BuD = DoubleSidedCrystalballFunction_BuD(x, par, bin_x, bin_y);
  //cout << DCB_BuD << " DCB_DuB\n";

  //auto signal = [this, bin_x, bin_y] (double *xx, double *param)->double { 
  //       return gausPDF(xx, param, bin_x, bin_y)+crystalballPDF(xx, param, bin_x, bin_y);
         //return crystalballPDF2(xx, param, bin_x, bin_y);
  //};
  //TF2 fsignal("fsig", signal, minx, maxx, miny, maxy, 87);
   //fsignal.SetParameters(par);
   //cout << fsignal.Integral(minx, maxx, miny, maxy)/bin_width/double(nevents) << " ######## intsig\n";
  //expo*=bin_width*double(nevents)*fcomb/int_expo;
  //double bin_width = (maxx-minx)/double(nbins);
  //cout << expo << " expo " << gaus1 << " gaus1 " << gaus2 << " CB\n";
  //cout <<int_expo<< "  " <<  expo << endl;
  //cout << CB << " " << gaus << "  " << expo << endl;
  //return hdata.Integral()*(fcomb*expo+(1.0-fcomb)*(f12*CB+gaus));
  //return bin_width*double(nevents)*((1.0-fcomb)*(f12*gausPDF(&m_rec, par)+(1.0-f12)*crystalballPDF(&m_rec, &par[0]))
                                       //+fcomb*expo);
  //return (expo+gaus1+gaus2+CBs+gausD);
  //return (cheb+gaus1+gaus2);
  //return (cheb+gaus1+gaus2+CBs)*bin_widthy;
  return (cheb+gaus1+gaus2+CBs+gausB0D+gausBuD+DCB_B0D+DCB_BuD)*bin_widthy;
  //return (expo+gaus1+gaus2+CBs+gausD+DCB);
  //return bin_width*nevents*((1.0-fcomb)*(f12*gausPDF(&m_rec, par)+(1.0-f12)*crystalballPDF(&m_rec, par))+fcomb*exponentPDF(&m_rec, par));
}
double Dmass::calcChi2(const double *param ){

    //for (int i=0; i<8; i++) cout << param[i] << "  ";
    //cout << endl;
    double chi2 = 0.;
    double loglike = 0.;
    double bin_width = (maxx-minx)/double(nbins);
    double bin_widthy = (maxy-miny)/double(nbins);
    auto nevents = vect_Dmass.size();
      for ( const auto& ivect: vect_Dmass){
	double mass[2];
        mass[0] = get<0>(ivect);
	mass[1] = get<1>(ivect);
        double bincont = Dmass_func(mass, param)/(bin_width*bin_widthy*double(nevents));
	if (bincont>1.0) cout << bincont << " <<<<<<<<<<<<<<<<<<<<<<<\n";
	if (bincont<=0.0) continue;
        //cout<< bin_width << "  " << double(nevents) << "  " << Dmass_func(mass, param) << endl;
        //double diff = bincont-entries;
      
      	double f12 = param[4];
        double fcomb = abs(param[5]);
        double fCB = abs(param[8]);
        double fDCB = abs(param[18]);
	double fBu = abs(param[23]);
	//err = sqrt(err*err+errMC2);
        //if (err==0) continue;
        //chi2+= diff*diff/err/err;
	//cout << chi2 << "  " << err << endl;
        //double fact_sum = 0.0;
        loglike-=log(bincont); 
	//cout << loglike << " SKlsnclancdjscnlskjncnjskn vnfjnvlnfj\n";
        //cout << nbins << " bins " << err << " err " << diff << " diff\n";
        //cout << bincont << " fval " << entries << " entries " << chi2 << " chi2++ \n";
//	cout << bincont << "  " << loglike << endl;
      }
    
        loglike*=2.0;//double(nevents);

	        for (int i=0; i<ncontr; i++){
                       for (int j=0; j<nprm; j++){
			          //cout << xx[i][j] << " tutu " << dxx[i][j] << endl; 
                                  double par2 = param[28+i*nprm+j];
				  double dXX = dxx[i][j];
                                  if(dXX>0){
					  double XX = xx[i][j];
                                          loglike+=((par2-XX)*(par2-XX))/(dXX*dXX);
                                  }
                       }
                }
		
		for (int i=0; i<nsigBu+nBs+nB0D+nBuD; i++){
			double par2 = param[i];
			          //cout << yy[i] << " tutu yy " << dyy[i] << endl; 
			if (dyy[i]>0) loglike+=((par2-yy[i])*(par2-yy[i]))/(dyy[i]*dyy[i]);
		}
		

    return loglike;
}
double Dmass::gaussian_error(double mean, double error, const double *param, int ipar){
	double p = param[ipar];
	error*=1.0;
	double chi2 = (p-mean)*(p-mean)/(error*error);
	return chi2;	
}


void Dmass::Draw(const double *param, const char *name, double CovMtx[nvar*nvar]){
    TCanvas *cDmass = new TCanvas(name, "Dmass", 1000, 500);
    TH2D *hdrawdata = dynamic_cast<TH2D*>(hdata.Clone(Form("%shdraw", name)));

    hdrawdata->SetStats(kFALSE);
    //hdrawdata->GetYaxis()->SetRangeUser(1e-3, 200000);

    TString s;
    if (sign) s = " M(K^{+}#pi^{-}#pi^{-})";
    else s = " M(K^{-}#pi^{+}#pi^{+})";
    hdrawdata->SetTitle(s);
    hdrawdata->SetMarkerStyle(kFullCircle);
    cDmass->cd();
    hdrawdata->DrawClone("ep");
    //pad1->SetLogz();
    //gPad->SetLogz();
    //TH2D hdummy;
    Dmass_autofit DmassDraw(vect_Dmass, hcorr, minx, maxx, miny, maxy);
    double xxx[2] = {1870., 5000.};
    double par[87];
    for (int ipar = 0; ipar<87; ipar++){
	    par[ipar] = param[ipar];
    }


    auto lambda_draw = [this](double *x, double *ppar)->double{ return Dmass_func(x, ppar);};

   // TF2 *fdraw = new TF2(Form("%sfdraw", name), DmassDraw, minx, maxx, miny, maxy, 87);
    TF2 *fdraw = new TF2(Form("%sfdraw", name), lambda_draw, minx, maxx, miny, maxy, 87);
    fdraw->SetNpx(1000);
    cout << fdraw << " fdraw \n\n";
    fdraw->SetLineColor(kRed);
    fdraw->SetParameters(param);
    cDmass->cd();
    fdraw->DrawClone("surf same");
    cDmass->SaveAs("fit2D.pdf");
    cDmass->SaveAs("fit2D.C");
    cDmass->SaveAs("fit2D.root");
    cout << "==============================\n " << Dmass_func(xxx, par ) << endl;
    cout << " eval fdraw "<< fdraw->Eval(1870., 5000.) << endl; 
    cout << fdraw->Integral(minx, maxx, miny, maxy)/bin_width/bin_widthy/double(nevents)<< "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< int " << endl;
    double pa[ncontr-1];
    for (int ipar=0; ipar<ncontr-1; ipar++){
	    pa[ipar] = fdraw->GetParameter(82+ipar);
    }
    double frac[ncontr];
    frac[0] = 1-pa[0];
    frac[1] = pa[0]*(1.0-pa[1]);
    frac[2] = pa[0]*pa[1]*(1.0-pa[2]);
    frac[3] = pa[0]*pa[1]*pa[2]*(1.0-pa[3]);
    frac[4] = pa[0]*pa[1]*pa[2]*pa[3]*(1.0-pa[4]);
    frac[5] = pa[0]*pa[1]*pa[2]*pa[3]*pa[4];

       double vect[ncontr][ncontr-1] = { -1.0, 0.0, 0.0, 0.0, 0.0,
                                   1.0-abs(pa[1]), -abs(pa[0]), 0.0, 0.0, 0.0,
                                   abs(pa[1])*(1.0-abs(pa[2])), abs(pa[0])*(1.0-abs(pa[2])), -abs(pa[0])*abs(pa[1]), 0.0, 0.0,
                                   abs(pa[1])*abs(pa[2])*(1.0-abs(pa[3])), abs(pa[0])*abs(pa[2])*(1.0-abs(pa[3])), abs(pa[0])*abs(pa[1])*(1.0-abs(pa[3])), -abs(pa[0])*abs(pa[1])*abs(pa[2]), 0.0,
                                   abs(pa[1])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4])), abs(pa[0])*abs(pa[2])*abs(pa[3])*(1.0-abs(pa[4])), abs(pa[0])*abs(pa[1])*abs(pa[3])*(1.0-abs(pa[4])), abs(pa[0])*abs(pa[1])*abs(pa[2])*(1.0-abs(pa[4])), -abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3]),
                                   abs(pa[1])*abs(pa[2])*abs(pa[3])*abs(pa[4]), abs(pa[0])*abs(pa[2])*abs(pa[3])*abs(pa[4]), abs(pa[0])*abs(pa[1])*abs(pa[3])*abs(pa[4]), abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[4]), abs(pa[0])*abs(pa[1])*abs(pa[2])*abs(pa[3])};
   double error[ncontr];
   double sum[ncontr][ncontr];
   for (int ipar = 0; ipar<ncontr ; ipar++){
           for (int jpar=0; jpar<ncontr-1; jpar++){
                sum[ipar][jpar] = 0.0;
                for (int kpar=0; kpar<ncontr-1; kpar++){
                        sum[ipar][jpar] += vect[ipar][kpar]*CovMtx[(82+jpar)*nvar+82+kpar];
                        //cout << sum[ipar][jpar] << " sum  " << ipar << "  " << jpar << endl;
                }
           }
           double sum2 = 0.0;
           for (int mpar=0; mpar<ncontr-1; mpar++){
                   sum2+=sum[ipar][mpar]*vect[ipar][mpar];
           }
           error[ipar] = sqrt(sum2);
   }

   		TString names[ncontr] = {"sig", "Bu", "Bs", "B0D", "comb", "BuD"};
                for (int ic=0; ic<ncontr; ic++) cout << names[ic] << "  frac = " << frac[ic] << " +- "<< error[ic] <<  endl;


//    TH2D checkhist2(hdata);
    TH2D hpull(hdata);
    hpull.SetName("hpull");
    TH1D hpull1D("hpull1D", "", nbins, -8., 8.);
    for (int binx=1; binx<=hdata.GetNbinsX(); ++binx){
      for (int biny=1; biny<=hdata.GetNbinsY(); ++biny){
        double mass[2];
        //mass[0] = minx+(binx-0.5)*bin_width;
        mass[0] = hdata.GetXaxis()->GetBinCenter(binx);
        //mass[1] = miny+(biny-0.5)*bin_widthy;
        mass[1] = hdata.GetYaxis()->GetBinCenter(biny);
        //double estim = fdraw->Eval(mass[0], mass[1]);
        double estim = Dmass_func(mass, param);
	//checkhist2.SetBinContent(binx, biny, bincont);
	double entries = hdata.GetBinContent(binx, biny);
	double err_data = hdata.GetBinError(binx, biny);
        if (err_data==0) continue;
	//cout << entries << "  " << bincont << endl;
	hpull.SetBinContent(binx, biny, (entries-estim)/err_data);
	hpull1D.Fill((entries-estim)/err_data);
      }
    }
    //cout << checkhist2.Integral() << " checkhist2 integral " << hdata.Integral() << endl;


/*
    TH2D checkhist(hdata);
    TH2D hgaus2D(hdata);
    TH2D hCB2D(hdata);
    TH2D hBs2D(hdata);
    TH2D hDCB2D(hdata);
    TH2D hBu2D(hdata);
    TH2D hcomb2D(hdata);
    for (int binx=1; binx<=hdata.GetNbinsX(); ++binx){
      for (int biny=1; biny<=hdata.GetNbinsY(); ++biny){
        double mass[2];
        mass[0] = minx+(binx-0.5)*bin_width;
        mass[1] = miny+(biny-0.5)*bin_widthy;
	int bincorr_x, bincorr_y;
	GetCorrBin(bincorr_x, bincorr_y, mass[0], mass[1]);
        double bincont = Dmass_func(mass, param)*bin_widthy;
	checkhist.SetBinContent(binx, biny, bincont);
      }
    }
    cout << checkhist.Integral() << " checkhist integral " << hdata.Integral() << endl;
*/
    TCanvas *cpull = new TCanvas("cpull", "", 500, 500);
    hpull.DrawClone("surf");
    
    cpull->cd();
    cpull->SaveAs("pull2D.pdf");
    cpull->SaveAs("pull2D.C");
    cpull->SaveAs("pull2D.root");

    TCanvas *cpull1D = new TCanvas("cpull1D", "", 500, 500);
    hpull1D.DrawClone();
    cpull1D->SaveAs("pull1D.pdf");
    cpull1D->SaveAs("pull1D.C");
}

Dmass_autofit::Dmass_autofit(std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_autofit::operator()(double *x, double *param){
    return Dmass_func(x, param);
}
Dmass_genExpo::Dmass_genExpo(std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genExpo::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return exponentPDF(x, param, bin_x, bin_y);
}
Dmass_genCheb::Dmass_genCheb(std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genCheb::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return chebyshevPDF(x, param, bin_x, bin_y);
}

Dmass_genGaus::Dmass_genGaus(std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genGaus::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return gausPDF(x, param, bin_x, bin_y);
}


Dmass_genCB2::Dmass_genCB2(std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genCB2::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return crystalballPDF2(x, param, bin_x, bin_y);
}
Dmass_genD0B::Dmass_genD0B(std::vector<std::pair<double, double>> vect_Dmass_, TH2D** hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genD0B::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return DoubleSidedCrystalballFunction_B0D(x, param, bin_x, bin_y)+gausPDF_B0D(x, param, bin_x, bin_y);
}
Dmass_genDuB::Dmass_genDuB(std::vector<std::pair<double, double>> vect_Dmass_, TH2D** hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genDuB::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return DoubleSidedCrystalballFunction_BuD(x, param, bin_x, bin_y)+gausPDF_BuD(x, param, bin_x, bin_y);
}
Dmass_genGaus2::Dmass_genGaus2(std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_): Dmass(vect_Dmass_, hcorr_, minx_, maxx_, miny_, maxy_){}
double Dmass_genGaus2::operator()(double *x, double *param){
    int bin_x, bin_y;
    GetCorrBin(bin_x, bin_y, x[0], x[1]);    
    return crystalballPDF(x, param, bin_x, bin_y);
}

void Dmass::normalize_histogram(TH2D &hist){
   int nbins_corr = hist.GetNbinsX();
   double nbinsx = double(nbins_corr);
   for (int binx=1; binx<=nbins_corr; binx++){
        double sum = 0.0;
        for (int biny=1; biny<=nbins; biny++){
          sum+=double(hist.GetBinContent(binx, biny));
        }
        if (sum==0) {
                nbinsx -= 1.;
                continue;
        }
        for (int biny=1; biny<=nbins; biny++){
          hist.SetBinContent(binx, biny, double(hist.GetBinContent(binx, biny))/sum);
          hist.SetBinError(binx, biny, double(hist.GetBinError(binx, biny))/sum);
        }
   }
   for (int binx=1; binx<=nbins_corr; binx++){
     double sum = 0.0;
     for (int biny=1; biny<=nbins; biny++){
          sum+=double(hist.GetBinContent(binx, biny));
     }
     if (sum==0) continue;
     for (int biny=1; biny<=nbins; biny++){
       if (nbinsx!=0) hist.SetBinContent(binx, biny, double(hist.GetBinContent(binx, biny))*double(nbins_corr)/nbinsx);
       if (nbinsx!=0) hist.SetBinError(binx, biny, double(hist.GetBinError(binx, biny))*double(nbins_corr)/nbinsx);
     }
  }
}
  

void Dmass::generate(const double *param, int N){
    //gRandom->SetSeed(0);

    double f12 = param[6];
    double fcomb = param[7];
    Dmass_genExpo GenExpo(vect_Dmass, hcorr, minx, maxx, miny, maxy);
    TF2 f1_expo("f1_expo", GenExpo, minx, maxx, miny, maxy, 87);
    f1_expo.SetParameters(param);
    //cout << f1_expo.Integral(minx, maxx) << " exp integral\n";

    Dmass_genGaus GenGaus(vect_Dmass, hcorr, minx, maxx, miny, maxy); 
    TF2 f1_gaus("f1_gaus", GenGaus, minx, maxx, miny, maxy, 87);
    f1_gaus.SetParameters(param); 
    //cout << f1_gaus.Integral(minx, maxx) << " gaus integral\n";
    Dmass_genGaus2 GenGaus2(vect_Dmass, hcorr, minx, maxx, miny, maxy); 
    TF2 f1_gaus2("f1_gaus2", GenGaus2, minx, maxx, miny, maxy, 87);
    f1_gaus2.SetParameters(param); 
    //cout << f1_gaus2.Integral(minx, maxx) << " gaus2 integral\n";
    for (int i=0; i<=std::round(fcomb*double(N)); ++i){
       double xx, yy;
       f1_expo.GetRandom2(xx, yy);	    
       hdata.Fill(xx, yy);
    }
    int exp_int = hdata.Integral();
    cout << hdata.Integral() << " integral data \n";
    for (int i=0; i<=std::round((1.0-fcomb)*f12*double(N)); ++i){
       double xx, yy;
       f1_gaus.GetRandom2(xx, yy);	    
       hdata.Fill(xx, yy);
    }
    int gaus1_int = hdata.Integral()-exp_int;
    cout << hdata.Integral() << " integral data \n";
    for (int i=0; i<=std::round((1.0-fcomb)*(1.0-f12)*double(N)); ++i){
       double xx, yy;
       f1_gaus2.GetRandom2(xx, yy);	    
       hdata.Fill(xx, yy);
    }
    int gaus2_int = hdata.Integral()-exp_int-gaus1_int;
    cout << hdata.Integral() << " integral data \n";

   //cout << "fcomb: " << double(exp_int)/double(exp_int+gaus1_int+gaus2_int) <<  " f12 " << double(gaus1_int)/double(gaus2_int)<< endl;
   
}

}
