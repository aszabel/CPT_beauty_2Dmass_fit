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

#include <iostream>

namespace cpt_b0_analysis{

double MMmin = 1800.;
double MMmax = 1950.;
double Mcorrmin = 3500.;
double Mcorrmax= 6000.;


double gaus_pdf (double *x, double *par){
  int nbins = 20;
  double bin_width = 12./int(sqrt(nbins));
  double sigma = par[0];
  double mean = par[1];
  return double(nbins)*bin_width*ROOT::Math::gaussian_pdf(x[0], sigma, mean);
}

Dmass::Dmass(TH2D &hdata_, TH2D **hcorr_, bool sign_) : hdata(hdata_), hcorr(hcorr_), sign(sign_){
   nbins = hdata.GetNbinsX();
   nevents = hdata.Integral();
   bin_width = (MMmax-MMmin)/double(nbins);
   bin_widthy = (Mcorrmax-Mcorrmin)/double(nbins);
   cout << hdata.GetYaxis()->GetBinLowEdge(15) <<  " maxbin " << endl;
}

Dmass::~Dmass(){}

double Dmass::operator()(const double *par){
  //cout << " OP \n";
  double chi2 = calcChi2(par); 
  //cout << chi2 << " chi2_result " << endl;
  return chi2;

}
double Dmass::gausPDF(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = par[0];
  double mean = par[1];
  double f12 = abs(par[4]);
  double fcomb = abs(par[5]);
  double fCB = abs(par[8]);
  double fD = abs(par[16]);
  double fDCB = abs(par[18]);
  double fBu = abs(par[23]);
  //double nevents = par[24];
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(MMmax, sigma, mean)-ROOT::Math::normal_cdf(MMmin, sigma, mean);
  gaus/=int_gaus;
  TH2D hist(*hcorr[0]);
  hist.Scale(1.0-fBu);
  hist.Add(hcorr[1], fBu);
  //hist.Scale(1./hist.Integral());
  normalize_histogram(hist);
  auto global_bin = hist.FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hist.GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hist.GetBinContent(bin_x, bin_y);
  
  return bin_width*double(nevents)*abs(1.0-fcomb)*abs(1.0-fD-fDCB-fCB)*f12*gaus*cont;
}
double Dmass::gausPDF_DspDm(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = par[15];
  double mean = par[1];
  double f12 = par[4];
  double fD = par[16];
  double fcomb = abs(par[5]);
  double gaus = ROOT::Math::gaussian_pdf(m_rec, sigma, mean);
  double int_gaus =  ROOT::Math::normal_cdf(MMmax, sigma, mean)-ROOT::Math::normal_cdf(MMmin, sigma, mean);
  gaus/=int_gaus;
  auto global_bin = hcorr[1]->FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hcorr[1]->GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hcorr[1]->GetBinContent(bin_x, bin_y);

  return bin_width*double(nevents)*(1.0-fcomb)*fD*gaus*cont;
}
double Dmass::crystalballPDF(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = par[0];
  double mean = par[1];
  double sigma2 = par[3];
  double f12 = abs(par[4]);
  double fcomb = abs(par[5]);
  double alpha = par[6];
  double n = par[7];
  double fCB = abs(par[8]);
  double fD = abs(par[16]);
  double fDCB = abs(par[18]);
  double fBu = abs(par[23]);
  //double nevents = par[24];
  //double CB = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean);
  //double int_CB =  ROOT::Math::normal_cdf(MMmax, sigma2, mean)-ROOT::Math::normal_cdf(MMmin, sigma2, mean);
  double CB = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma2, mean);
  double int_CB = TMath::Abs(ROOT::Math::crystalball_integral(MMmin, alpha, n, sigma2, mean)-ROOT::Math::crystalball_integral(MMmax, alpha, n, sigma2, mean));
  CB/=int_CB;

  TH2D hist(*hcorr[0]);
  hist.Scale(1.0-fBu);
  hist.Add(hcorr[1], fBu);  
  //hist.Scale(1./hist.Integral());
  normalize_histogram(hist);
  auto global_bin = hist.FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hist.GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hist.GetBinContent(bin_x, bin_y);
  return bin_width*double(nevents)*abs(1.0-fcomb)*abs(1.0-fD-fCB-fDCB)*abs(1.0-f12)*CB*cont;
}
double Dmass::DoubleSidedCrystalballFunction(double *x, const double *par)
{
  double alpha = par[20];
  double n     = par[21];
  double mean  = par[1];
  double sigma = par[19];
  double fcomb = par[5];
  double f12 = par[4];
  double fCB = par[8];
  double fDCB = abs(par[18]);
  double alpha_h = par[22];
  double m_rec=x[0];
  double m_corr = x[1];
  double result;

   if (m_rec<mean){
       result = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma, mean);
   }
   else{
       result = ROOT::Math::crystalball_function(2.*mean-m_rec, alpha_h, n, sigma, mean);
   }
      double intDCB = TMath::Abs(-ROOT::Math::crystalball_integral(MMmin, alpha, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha, n, sigma, mean)) + TMath::Abs(-ROOT::Math::crystalball_integral(2.*mean-MMmax, alpha_h, n, sigma, mean)+ROOT::Math::crystalball_integral(mean, alpha_h, n, sigma, mean));
  
  auto global_bin = hcorr[3]->FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hcorr[3]->GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hcorr[3]->GetBinContent(bin_x, bin_y);
   return bin_width*double(nevents)*(1.0-fcomb)*fDCB*result/intDCB*cont;
}

double Dmass::crystalballPDF2(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = par[0];
  double mean = par[9];
  double sigma2 = par[10];
  double f12 = par[4];
  double fcomb = abs(par[5]);
  double alpha = par[11];
  double n = par[12];
  double fCB = abs(par[8]);
  double fDCB = abs(par[18]);
  //double CB = ROOT::Math::gaussian_pdf(m_rec, sigma2, mean);
  //double int_CB =  ROOT::Math::normal_cdf(MMmax, sigma2, mean)-ROOT::Math::normal_cdf(MMmin, sigma2, mean);
  double CB = ROOT::Math::crystalball_function(m_rec, alpha, n, sigma2, mean);
  double int_CB = TMath::Abs(ROOT::Math::crystalball_integral(MMmin, alpha, n, sigma2, mean)-ROOT::Math::crystalball_integral(MMmax, alpha, n, sigma2, mean));
  CB/=int_CB;
  //return bin_width*double(nevents)*(1.0-fcomb)*(1.0-f12)*abs(1.0-fCB)*abs(1.0-fDCB)*CB;
  
  auto global_bin = hcorr[2]->FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hcorr[2]->GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hcorr[2]->GetBinContent(bin_x, bin_y);
  return bin_width*double(nevents)*(1.0-fcomb)*fCB*CB*cont;
}
double Dmass::exponentPDF(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double exp_slope = par[2];
  double fcomb = abs(par[5]);
  //double expo = TMath::Exp(-exp_slope*m_rec);
  //double int_expo = 1.0/exp_slope *(TMath::Exp(-exp_slope*MMmin)-TMath::Exp(-exp_slope*MMmax));
  double expo = ROOT::Math::exponential_pdf(m_rec, exp_slope);
  double int_expo =  ROOT::Math::exponential_cdf(MMmax, exp_slope)-ROOT::Math::exponential_cdf(MMmin, exp_slope);
  expo /= int_expo;
 
  auto global_bin = hcorr[4]->FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hcorr[4]->GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hcorr[4]->GetBinContent(bin_x, bin_y);
  return bin_width*double(nevents)*fcomb*expo*cont;
}
double Dmass::chebyshevPDF(double *x, const double *par){
  double m_rec = x[0];
  double m_corr = x[1];
  double fcomb = abs(par[5]);
  double a1 = par[13];
  double a2 = par[14];
  //double nevents = par[24];
  double cheb = 1.0 + a1*m_rec + a2*(2.0*m_rec*m_rec-1.0);
  double int_cheb = (1.0-a2)*MMmax + 0.5*a1*MMmax*MMmax+2./3.*a2*MMmax*MMmax*MMmax- (1.0-a2)*MMmin-0.5*a1*MMmin*MMmin-2./3.*a2*MMmin*MMmin*MMmin;
  cheb/=int_cheb;
  
  auto global_bin = hcorr[4]->FindBin(m_rec, m_corr);
  int bin_x, bin_y, bin_z;
  hcorr[4]->GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
  double cont  = hcorr[4]->GetBinContent(bin_x, bin_y);
  return bin_width*double(nevents)*fcomb*cheb*cont;
}
double Dmass::Dmass_func(double *x, const double *par){

  double m_rec = x[0];
  double m_corr = x[1];
  double sigma = par[0];
  double mean  = par[1];
  double exp_slope = par[2];
  double s12 = par[3];
  double f12 = par[4];
  double fcomb = par[5];
  double alpha = par[6];
  double n = par[7];

                              
  double nevents = double(hdata.Integral());
  if (nevents==0) nevents = 1000000.;

  double expo = exponentPDF(x, par);
  //cout << expo << "expo\n";
  double gaus1 = gausPDF(x, par);
  //cout << gaus1 << "gaus1\n";
  double gausD = gausPDF_DspDm(x, par);
  //cout << gausD << "gausD\n";
  double gaus2 =  crystalballPDF(x, par);
  //cout << gaus2 << "gaus2\n";
  double CBs = crystalballPDF2(x, par);
  //cout << CBs << "CBs\n";
  double cheb = chebyshevPDF(x, par);
  //cout << cheb << "cheb\n";
  double DCB = DoubleSidedCrystalballFunction(x, par);
  //cout << DCB << "DCB\n";

  //expo*=bin_width*double(nevents)*fcomb/int_expo;
  //double bin_width = (MMmax-MMmin)/double(nbins);
  //cout << expo << " expo " << gaus1 << " gaus1 " << gaus2 << " CB\n";
  //cout <<int_expo<< "  " <<  expo << endl;
  //cout << CB << " " << gaus << "  " << expo << endl;
  //return hdata.Integral()*(fcomb*expo+(1.0-fcomb)*(f12*CB+gaus));
  //return bin_width*double(nevents)*((1.0-fcomb)*(f12*gausPDF(&m_rec, par)+(1.0-f12)*crystalballPDF(&m_rec, &par[0]))
                                       //+fcomb*expo);
  //return (expo+gaus1+gaus2+CBs+gausD);
  //return (cheb+gaus1+gaus2);
  return (cheb+gaus1+gaus2+CBs+gausD+DCB);
  //return (expo+gaus1+gaus2+CBs+gausD+DCB);
  //return bin_width*nevents*((1.0-fcomb)*(f12*gausPDF(&m_rec, par)+(1.0-f12)*crystalballPDF(&m_rec, par))+fcomb*exponentPDF(&m_rec, par));
}
double Dmass::calcChi2(const double *param ){

    //for (int i=0; i<8; i++) cout << param[i] << "  ";
    //cout << endl;
    double chi2 = 0.;
    double loglike = 0.;
    double bin_width = (MMmax-MMmin)/double(nbins);
    double bin_widthy = (Mcorrmax-Mcorrmin)/double(nbins);
    for (int binx=1; binx<=hdata.GetNbinsX(); ++binx){
      for (int biny=1; biny<=hdata.GetNbinsY(); ++biny){
	double mass[2];
        mass[0] = MMmin+(binx-0.5)*bin_width;
	mass[1] = Mcorrmin+(biny-0.5)*bin_widthy;
        double bincont = Dmass_func(mass, param);
        double entries = hdata.GetBinContent(binx, biny);
        double diff = bincont-entries;
        double err = hdata.GetBinError(binx, biny);
      
      	double f12 = param[4];
        double fcomb = abs(param[5]);
        double fCB = abs(param[8]);
        double fDCB = abs(param[18]);
	double fBu = abs(param[23]);
	//double nevents = param[24];
        TH2D htemp(*hcorr[0]);
	auto global_bin = htemp.FindBin(mass[0], mass[1]);
        int bin_x, bin_y, bin_z;
        htemp.GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
	double errMC = htemp.GetBinError(bin_x, bin_y);
	TH2D htempBu(*hcorr[1]);
	global_bin = htempBu.FindBin(mass[0], mass[1]);
        htempBu.GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
        double errBu = htempBu.GetBinError(bin_x, bin_y);
	errMC = sqrt((1.0-fBu)*(1.0-fBu)*errMC*errMC+fBu*fBu*errBu*errBu);
	htemp.Scale(1.0-fBu);
	htemp.Add(&htempBu, fBu);
        double errgaus, errCB1;
	if (htemp.GetBinContent(bin_x, bin_y)==0.0) {
		errgaus = 0.0;
		errCB1 = 0.0;
	}else{
		errgaus = gausPDF(mass, param)/htemp.GetBinContent(bin_x, bin_y)*errMC;
		errCB1 = crystalballPDF(mass, param)/htemp.GetBinContent(bin_x, bin_y)*errMC;
	}
	errMC = sqrt(errgaus*errgaus+errCB1*errCB1);
        //htemp.Scale(1./htemp.Integral());
	TH2D htempCB(*hcorr[2]);
	double errCB = htempCB.GetBinError(bin_x, bin_y);
	if (htempCB.GetBinContent(bin_x, bin_y)==0.0) errCB = 0.0;
	else errCB *= crystalballPDF2(mass, param)/htempCB.GetBinContent(bin_x, bin_y);

	TH2D htempDCB(*hcorr[3]);
	double errDCB = htempDCB.GetBinError(bin_x, bin_y);
	if (htempDCB.GetBinContent(bin_x, bin_y)==0) errDCB=0.0;
	else errDCB *= DoubleSidedCrystalballFunction(mass, param)/htempDCB.GetBinContent(bin_x, bin_y);

	TH2D htempComb(*hcorr[4]);
	double errComb = htempComb.GetBinError(bin_x, bin_y);
	if (htempComb.GetBinContent(bin_x, bin_y)==0.0) errComb = 0.0;
	else errComb *= chebyshevPDF(mass, param)/htempComb.GetBinContent(bin_x, bin_y);
	double errMC2 = (errMC*errMC+errCB*errCB+errDCB*errDCB+errComb*errComb);

	err = sqrt(err*err+errMC2);
        if (err==0) continue;
        chi2+= diff*diff/err/err;
	//cout << chi2 << "  " << err << endl;
        double fact_sum = 0.0;
        for (int i=1; i<=entries; ++i){
          fact_sum+=log(i);
        }
        loglike+= bincont - entries*log(bincont)+fact_sum;
        //cout << nbins << " bins " << err << " err " << diff << " diff\n";
        //cout << bincont << " fval " << entries << " entries " << chi2 << " chi2++ \n";
      }
    }
        loglike*=2.0;

    double ext_chi2 = gaussian_error(1879., 1.2, param, 9);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(13.5, 0.72, param, 10);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(0.370, 0.066, param, 11);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(12.43, 21., param, 12);
    chi2 += ext_chi2;
    loglike += ext_chi2;

    chi2+=gaussian_error(7.41, 0.20, param, 15);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(1870.5, 0.23, param, 17);
    chi2 += ext_chi2;

    loglike += ext_chi2;
    ext_chi2 = gaussian_error(7.204, 0.134, param, 19);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(1.8423, 0.092, param, 20);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(2.183, 0.332, param, 21);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    ext_chi2 = gaussian_error(1.8423, 0.092, param, 22);
    chi2 += ext_chi2;
    loglike += ext_chi2;
    //return chi2/double(nbins);
    return chi2;
    //return loglike;
}
double Dmass::gaussian_error(double mean, double error, const double *param, int ipar){
	double p = param[ipar];
	double chi2 = (p-mean)*(p-mean)/(error*error);
	return chi2;	
}

void Dmass::Draw(const double *param, const char *name){    
    TCanvas *cDmass = new TCanvas(name, "Dmass", 1000, 500);
    //cDmass->Divide(2, 1);
    //cDmass->cd(1);
    //TPad *pad1 = new TPad(Form("%spad1", name), "", 0.0, 0.3, 1.0, 1.0);
    //pad1->Draw();
    //pad1->cd(); 
    TH2D *hdrawdata = dynamic_cast<TH2D*>(hdata.Clone(Form("%shdraw", name)));
    //hdrawdata->SetMarkerStyle(kFullCircle);
    //hdrawdata->SetMarkerSize(0.5);

    hdrawdata->SetStats(kFALSE);
    hdrawdata->GetYaxis()->SetRangeUser(10, 200000);

    TString s;
    if (sign) s = " M(K^{+}#pi^{-}#pi^{-})";
    else s = " M(K^{-}#pi^{+}#pi^{+})";
    hdrawdata->SetTitle(s);
    hdrawdata->SetMarkerStyle(kFullCircle);
    hdrawdata->DrawClone("ep");
    //pad1->SetLogz();
    //gPad->SetLogz();

    //TH2D hdummy;
    Dmass_autofit DmassDraw(hdata, hcorr);
    TF2 *fdraw = new TF2(Form("%sfdraw", name), DmassDraw, MMmin, MMmax, Mcorrmin, Mcorrmax, 25);
    fdraw->SetLineColor(kRed);
    fdraw->SetParameters(param);

    TH2D checkhist2(hdata);
    TH2D hpull(hdata);
    TH1D hpull1D("hpull1D", "", nbins, -8., 8.);
    for (int binx=1; binx<=hdata.GetNbinsX(); ++binx){
      for (int biny=1; biny<=hdata.GetNbinsY(); ++biny){
        double mass[2];
        mass[0] = MMmin+(binx-0.5)*bin_width;
        mass[1] = Mcorrmin+(biny-0.5)*bin_widthy;
        double bincont = fdraw->Eval(mass[0], mass[1]);
	checkhist2.SetBinContent(binx, biny, bincont);
	double entries = hdata.GetBinContent(binx, biny);
	double err_data = hdata.GetBinError(binx, biny);
	double f12 = param[4];
        double fcomb = abs(param[5]);
        double fCB = abs(param[8]);
        double fDCB = abs(param[18]);
        double fBu = abs(param[23]);
        //double nevents = param[24];
	TH2D htemp(*hcorr[0]);
        auto global_bin = htemp.FindBin(mass[0], mass[1]);
        int bin_x, bin_y, bin_z;
        htemp.GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
        double errMC = htemp.GetBinError(bin_x, bin_y);
        TH2D htempBu(*hcorr[1]);
        global_bin = htempBu.FindBin(mass[0], mass[1]);
        htempBu.GetBinXYZ(global_bin, bin_x, bin_y, bin_z);
        double errBu = htempBu.GetBinError(bin_x, bin_y);
        errMC = sqrt((1.0-fBu)*(1.0-fBu)*errMC*errMC+fBu*fBu*errBu*errBu);
        htemp.Scale(1.0-fBu);
        htemp.Add(&htempBu, fBu);
	
	double errgaus, errCB1;
        if (htemp.GetBinContent(bin_x, bin_y)==0.0) {
                errgaus = 0.0;
                errCB1 = 0.0;
        }else{
                errgaus = gausPDF(mass, param)/htemp.GetBinContent(bin_x, bin_y)*errMC;
                errCB1 = crystalballPDF(mass, param)/htemp.GetBinContent(bin_x, bin_y)*errMC;
        }
        errMC = sqrt(errgaus*errgaus+errCB1*errCB1);
        //htemp.Scale(1./htemp.Integral());
        TH2D htempCB(*hcorr[2]);
        double errCB = htempCB.GetBinError(bin_x, bin_y);
        if (htempCB.GetBinContent(bin_x, bin_y)==0.0) errCB = 0.0;
        else errCB *= crystalballPDF2(mass, param)/htempCB.GetBinContent(bin_x, bin_y);

        TH2D htempDCB(*hcorr[3]);
        double errDCB = htempDCB.GetBinError(bin_x, bin_y);
        if (htempDCB.GetBinContent(bin_x, bin_y)==0) errDCB=0.0;
        else errDCB *= DoubleSidedCrystalballFunction(mass, param)/htempDCB.GetBinContent(bin_x, bin_y);

        TH2D htempComb(*hcorr[4]);
        double errComb = htempComb.GetBinError(bin_x, bin_y);
        if (htempComb.GetBinContent(bin_x, bin_y)==0.0) errComb = 0.0; 
        else errComb *= chebyshevPDF(mass, param)/htempComb.GetBinContent(bin_x, bin_y);
        double errMC2 = (errMC*errMC+errCB*errCB+errDCB*errDCB+errComb*errComb);


        err_data = sqrt(err_data*err_data+errMC2);
        if (err_data==0) continue;
	//cout << entries << "  " << bincont << endl;
	hpull.SetBinContent(binx, biny, (entries-bincont)/err_data);
	hpull1D.Fill((entries-bincont)/err_data);
	//if (abs((entries-bincont)/err_data)>4.0) {
	//	cout << binx << "  " << biny << "  " << (entries-bincont)/err_data << "  " << entries-bincont << "  " << err_data << "  " << sqrt(errMC2) << "  " << errMC << "  " << errCB << "  " << errDCB << "  " << errComb << endl;
 	//}
      }
    }
    cout << checkhist2.Integral() << " checkhist2 integral " << hdata.Integral() << endl;

    cout << " eval fdraw "<< fdraw->Eval(1870.) << endl; 
    fdraw->DrawClone("same surf");

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
        mass[0] = MMmin+(binx-0.5)*bin_width;
        mass[1] = Mcorrmin+(biny-0.5)*bin_widthy;
        double bincont = Dmass_func(mass, param);
        checkhist.SetBinContent(binx, biny, bincont);
        bincont = gausPDF(mass, param);
        hgaus2D.SetBinContent(binx, biny, bincont);
        bincont = crystalballPDF(mass, param);
        hCB2D.SetBinContent(binx, biny, bincont);
        bincont = crystalballPDF2(mass, param);
        hBs2D.SetBinContent(binx, biny, bincont);
        bincont = DoubleSidedCrystalballFunction(mass, param);
        hDCB2D.SetBinContent(binx, biny, bincont);
        bincont = gausPDF_DspDm(mass, param);
        hBu2D.SetBinContent(binx, biny, bincont);
        bincont = chebyshevPDF(mass, param);
        hcomb2D.SetBinContent(binx, biny, bincont);
      }
    }

    cout << checkhist.Integral() << " checkhist integral " << hdata.Integral() << endl;

    TCanvas *cpull = new TCanvas("cpull", "", 500, 500);
    hpull.DrawClone("surf");
    TCanvas *cpull1D = new TCanvas("cpull1D", "", 500, 500);
    hpull1D.DrawClone();


    TH1D *hfunc1D  = (TH1D*)checkhist.ProjectionX("hfunc1D");
    TH1D *hdata1D_DM = (TH1D*)hdata.ProjectionX("hdata1D_DM");
    TH1D *hgaus1D = (TH1D*)hgaus2D.ProjectionX("hgaus2D");
    TH1D *hCB1D = (TH1D*)hCB2D.ProjectionX("hCB1D");
    TH1D *hBs1D = (TH1D*)hBs2D.ProjectionX("hBs1D");
    TH1D *hDCB1D = (TH1D*)hDCB2D.ProjectionX("hDCB1D");
    TH1D *hBu1D = (TH1D*)hBu2D.ProjectionX("hBu1D_DM");
    TH1D *hcomb1D = (TH1D*)hcomb2D.ProjectionX("hcomb1D");

    TCanvas *c1D_DM = new TCanvas("c1D_DM", "", 500, 500);
    TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    hdata1D_DM->SetStats(kFALSE);
    hdata1D_DM->Draw("ep");
    hfunc1D->SetLineWidth(2);
    hfunc1D->Draw("same");
    hgaus1D->SetLineColor(kRed);
    hgaus1D->SetLineWidth(2);
    hgaus1D->Draw("same");
    hCB1D->SetLineColor(kBlack);
    hCB1D->SetLineWidth(2);
    hCB1D->Draw("same");
    hBs1D->SetLineColor(kMagenta);
    hBs1D->SetLineWidth(2);
    hBs1D->Draw("same");
    hDCB1D->SetLineColor(kGreen);
    hDCB1D->SetLineWidth(2);
    hDCB1D->Draw("same");
    hBu1D->SetLineColor(kOrange);
    hBu1D->SetLineWidth(2);
    hBu1D->Draw("same");
    hcomb1D->SetLineColor(kCyan);
    hcomb1D->SetLineWidth(2);
    hcomb1D->Draw("same");


    TLegend *leg = new TLegend(0.8, 0.5, 1.0, 1.0);
    leg->AddEntry(hdata1D_DM, "data", "pe");
    leg->AddEntry(hfunc1D, "fit", "l");
    leg->AddEntry(hgaus1D, "gaus", "l");
    leg->AddEntry(hCB1D, "CB", "l");
    leg->AddEntry(hDCB1D, "B0->DsmDp", "l");
    leg->AddEntry(hBu1D, "Bu->Ds0Dp", "l");
    leg->AddEntry( hBs1D, "Bs->Dsmunu", "l");
    leg->AddEntry(hcomb1D, "comb", "l");
    leg->Draw("same");

    hdata1D_DM->Sumw2();
    TH1D histpull1D(*hdata1D_DM);
    for (int bin=1; bin<=hdata.GetNbinsX(); bin++){
        double err = hdata1D_DM->GetBinError(bin);
        cout << err << endl;
        if (err==0) continue;
        double diff  = hdata1D_DM->GetBinContent(bin)-hfunc1D->GetBinContent(bin);
        cout << diff << " diff\n";
        histpull1D.SetBinContent(bin, diff/err);
    }
    c1D_DM->cd();
    TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
    pad2->Draw();
    pad2->cd();
    histpull1D.SetStats(kFALSE);
    histpull1D.SetFillColor(kBlue);
    histpull1D.GetYaxis()->SetLabelSize(0.1);
    histpull1D.GetXaxis()->SetLabelSize(0.1);
    histpull1D.DrawClone("hist");


/*    hExpo.SetLineColor(kGreen);
    hExpo.SetFillColor(kGreen);
    hExpo.DrawClone("hist same");
    hGaus.SetLineColor(kMagenta);
    hGaus.DrawClone("hist same");
    hCB.SetLineColor(kRed);
    hCB.DrawClone("hist same");
    TH2D *hsum = dynamic_cast<TH2D*> (hCB.Clone("hsum"));
    //hsum.Reset("ICESM");
    hsum->Add(&hGaus);
    hsum->Add(&hExpo);
    hsum->SetLineColor(kBlack);
    hsum->DrawClone("hist same");
     
    Dmass_genGaus gaus1(hdata, hcorr);
    TF1 fitgaus("fitgaus", gaus1, MMmin, MMmax, 25);
    fitgaus.SetParameters(param);
    fitgaus.SetLineColor(kBlue);
    //fitgaus.DrawClone("same");

   Dmass_genGausD gausD(hdata, hcorr);
    TF1 fitgausD("fitgausD", gausD, MMmin, MMmax, 25);
    fitgausD.SetParameters(param);
    fitgausD.SetLineColor(kOrange);
    //fitgausD.DrawClone("same");

    Dmass_genGaus2 gaus2(hdata, hcorr);
    TF1 fitgaus2("fitgaus2", gaus2, MMmin, MMmax, 25);
    fitgaus2.SetParameters(param);
    fitgaus2.SetLineColor(kBlack);
    //fitgaus2.DrawClone("same");

    Dmass_genDCB DCB(hdata, hcorr);
    TF2 fitDCB("fitDCB", DCB, MMmin, MMmax, Mcorrmin, Mcorrmax, 25);
    fitDCB.SetParameters(param);
    fitDCB.SetLineColor(kViolet);
    //fitDCB.DrawClone("same surf");
    //cout << fitDCB.Eval((MMmax-MMmin)/2., (Mcorrmax-Mcorrmin)/2.) << " DCB eval\n";
    cout << fitDCB.Eval(1869.49, (Mcorrmax-Mcorrmin)/2.) << " DCB eval\n";
    for (int i=0; i<5; i++){
    cout << hcorr[i]->GetBinContent(7, 1) << "  " << i << endl;
    }

    Dmass_genCB2 CB2(hdata, hcorr);
    TF1 fitCB2("fitCB2", CB2, MMmin, MMmax, 25);
    fitCB2.SetParameters(param);
    fitCB2.SetLineColor(kMagenta);
    //fitCB2.DrawClone("same");

    Dmass_genExpo Expofit(hdata, hcorr);
    TF1 fitexpo("fitexpo", Expofit, MMmin, MMmax, 25);
    fitexpo.SetParameters(param);
    fitexpo.SetLineColor(kGreen);
    fitexpo.SetFillColor(kGreen);
    //fitexpo.DrawClone("same");

    Dmass_genCheb Cheb(hdata, hcorr);
    TF1 fitcheb("fitcheb", Cheb, MMmin, MMmax, 25);
    fitcheb.SetParameters(param);
    fitcheb.SetLineColor(kGreen);
    fitcheb.SetFillColor(kGreen);
    //fitcheb.DrawClone("same");

  */
  
   /* cDmass->cd(1);
    TPad *pad2 = new TPad(Form("%spad2", name), "", 0.0, 0.0, 1.0, 0.3);
    pad2->Draw();
    pad2->cd();
*/

    //TH2D *hsum_pull = new TH2D (Form("%shsum", name), "", int(sqrt(nbins)), -6.0, 6.0);
    //TH2D hpull(Form("%shpull", name), "", nbins, hdrawdata->GetBinLowEdge(1), hdrawdata->GetBinLowEdge(nbins+1));
    /*double bin_width = (MMmax-MMmin)/double(nbins);
    for (int bin=1; bin<=nbins;++bin){
       double mass = MMmin+(bin-0.5)*bin_width;
       double diff = Dmass_func( &mass, param) -hdrawdata->GetBinContent(bin);
       double err = hdrawdata->GetBinError(bin);
       hpull.SetBinContent(bin, diff/err);
       hsum_pull->Fill(diff/err);
       //cout << this->Dmass_func( &mass, param) << "  " << fdraw->Eval(mass) << "  " << hdrawdata->GetBinContent(bin) 
       //     << "  " << hpull.GetBinContent(bin) << endl;
      // hpull.SetBinContent(bin, fdraw->Eval(mass) -hdrawdata->GetBinContent(bin));
       //cout << mass<<"  " << fdraw->Eval(mass)  << " pull " << hdrawdata->GetBinContent(bin) << endl;
    }
    hpull.SetStats(0);
    hpull.SetFillColor(kBlue);
    hpull.GetYaxis()->SetLabelSize(0.1);    
    hpull.DrawClone();



    cDmass->cd(2);

    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(10011);
    TF1 *fgaus = new TF1("fgaus", gaus_pdf, -4.0, 4.0, 2);
    fgaus->SetParameters(1.0, 0.0);
    fgaus->SetParNames("sigma_fit", "mean_fit");
    fgaus->SetParLimits(0, 0.0, 5.0);

    hsum_pull->DrawClone();
    hsum_pull->Fit("fgaus");
    cDmass->SaveAs("DmassFit.pdf");
    //delete hsum;
    //delete hdrawdata;
    //delete fdraw;
    */
}

Dmass_autofit::Dmass_autofit(TH2D& hdata_, TH2D **hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_autofit::operator()(double *x, double *param){
    return Dmass_func(x, param);
}
Dmass_genExpo::Dmass_genExpo(TH2D& hdata_, TH2D **hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genExpo::operator()(double *x, double *param){
    return exponentPDF(x, param);
}
Dmass_genCheb::Dmass_genCheb(TH2D& hdata_, TH2D **hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genCheb::operator()(double *x, double *param){
    return chebyshevPDF(x, param);
}

Dmass_genGaus::Dmass_genGaus(TH2D& hdata_, TH2D **hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genGaus::operator()(double *x, double *param){
    return gausPDF(x, param);
}

Dmass_genGausD::Dmass_genGausD(TH2D& hdata_, TH2D** hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genGausD::operator()(double *x, double *param){
    return gausPDF_DspDm(x, param);
}
Dmass_genCB2::Dmass_genCB2(TH2D& hdata_, TH2D **hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genCB2::operator()(double *x, double *param){
    return crystalballPDF2(x, param);
}
Dmass_genDCB::Dmass_genDCB(TH2D& hdata_, TH2D** hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genDCB::operator()(double *x, double *param){
    return DoubleSidedCrystalballFunction(x, param);
}
Dmass_genGaus2::Dmass_genGaus2(TH2D& hdata_, TH2D **hcorr_): Dmass(hdata_, hcorr_){}
double Dmass_genGaus2::operator()(double *x, double *param){
    return crystalballPDF(x, param);
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
    Dmass_genExpo GenExpo(hdata, hcorr);
    TF2 f1_expo("f1_expo", GenExpo, MMmin, MMmax, Mcorrmin, Mcorrmax, 25);
    f1_expo.SetParameters(param);
    //cout << f1_expo.Integral(MMmin, MMmax) << " exp integral\n";

    Dmass_genGaus GenGaus(hdata, hcorr); 
    TF2 f1_gaus("f1_gaus", GenGaus, MMmin, MMmax, Mcorrmin, Mcorrmax, 25);
    f1_gaus.SetParameters(param); 
    //cout << f1_gaus.Integral(MMmin, MMmax) << " gaus integral\n";
    Dmass_genGaus2 GenGaus2(hdata, hcorr); 
    TF2 f1_gaus2("f1_gaus2", GenGaus2, MMmin, MMmax, Mcorrmin, Mcorrmax, 25);
    f1_gaus2.SetParameters(param); 
    //cout << f1_gaus2.Integral(MMmin, MMmax) << " gaus2 integral\n";
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
