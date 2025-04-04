/**
 *  Adam Szabelski 13.06.2019 tool for generation and analysis of CPT conservation in semileptonic decays 
 *  A list of semileptonic decay time functions to be fitted and or generated
 *
 *
 *  @file pdfs_cpt_model.cpp
*/

#include "pdfs_cpt_model.h"

#include <TMath.h>
#include <iostream>
namespace cpt_b0_analysis
{


/// CPT sensitive ////////////////////////////////////////////
double pdfs::B_to_f(double *x, const double *par)
{
  return pdfs::calc_B_to_f(x, par, 1.0);
}
double pdfs::Bbar_to_fbar(double *x, const double *par)
{
  return pdfs::calc_B_to_f(x, par, -1.0);
}


double pdfs::calc_B_to_f(double *x, const double *par, double m) // (m=1 B_to_f) (m=-1 Bbar_to_fbar) 
{
  auto t = x[0];
  auto R = par[0];//*(-2.0)*Params::dM/Params::dGamma;
  auto I = par[1];
  auto Gamma = par[2];
  //auto Gamma = Params::Gamma;
  auto exp = TMath:: Exp(-Gamma*t);
  auto zeta =(R*R+I*I+1)*TMath::CosH(0.5*Params::dGamma*t)+(1-R*R-I*I)*TMath::Cos(Params::dM*t);
  auto Imz = 2.*I*TMath::Sin(Params::dM*t);
  auto Rez = 2.*R*TMath::SinH(0.5*Params::dGamma*t);

  return exp*(zeta-m*(Rez+Imz));///(pdfs::normalization_B_to_f(R, I, m)+pdfs::normalization_CP(R, I));
}


///////////// CP sensitive //////////////////////////////////////////////

double pdfs::calcB_to_fbar(double *x, const double *par, double m)// 
{
  auto t = x[0];
  auto R = par[0];//*(-2.0)*Params::dM/Params::dGamma;
  auto I = par[1];
  auto Gamma = par[2];
  //auto Gamma = Params::Gamma;
  auto exp = TMath:: Exp(-Gamma*t);

  auto zeta_CP = TMath::Sqrt((1.-R*R+I*I)*(1.-R*R+I*I)+4.*R*R*I*I)*(TMath::CosH(-0.5*Params::dGamma*t)-TMath::Cos(Params::dM*t));
  return exp*(zeta_CP);///(pdfs::normalization_B_to_f(R, I, m)+pdfs::normalization_CP(R, I));
}
double pdfs::Bbar_to_f(double *x, const double *par)// no omega^2 in the function. It is added in generetion.
{
  return pdfs::calcB_to_fbar(x, par, 1.0);
}
double pdfs::B_to_fbar(double *x, const double *par)// no omega^2 in the function. It is added in generetion.
{
  return calcB_to_fbar(x, par, -1.0);
}
//////////////////Fit functions//////////////////////////////////
//////////////////Fit functions//////////////////////////////////

double pdfs::fitf_unasym(double *x, double *par)
{
  auto eta=par[3];
  auto num = Params::N*((1.0+eta)*pdfs::B_to_f(x, par) + (1.0-eta)*pdfs::Bbar_to_f(x, par))-Params::Nbar*((1.0+eta)*pdfs::B_to_fbar(x, par) + (1.0-eta)*pdfs::Bbar_to_fbar(x, par));
  auto den = Params::N*((1.0+eta)*pdfs::B_to_f(x, par) + (1.0-eta)*pdfs::Bbar_to_f(x, par))+Params::Nbar*((1.0+eta)*pdfs::B_to_fbar(x, par) + (1.0-eta)*pdfs::Bbar_to_fbar(x, par));
  return num/den;
}
double pdfs::fitf_simult(double *x, double *par)
{
  auto eta=par[3];
  auto flavour = x[1];
  double sim;
  if(flavour <1.0) sim=Params::N*((1.0+eta)*pdfs::B_to_f(x, par) +(1.0-eta)*pdfs::Bbar_to_f(x, par));    //decays to fermion (flavour=0)
  else sim=Params::Nbar*((1.0+eta)*pdfs::B_to_fbar(x, par) + (1.0-eta)*pdfs::Bbar_to_fbar(x, par));//decays to antifermion (flavour=1)

  //return time_acceptance(x, &par[3])* sim; 
  return  sim;
}
double pdfs::time_acceptance(double *x, const double *par)
{
   auto t = x[0];
   //auto t_shift = Params::t_shift;
   //auto alpha = Params::alpha;
   //auto beta = Params::beta;
   double t_shift = par[0];
   double alpha = par[1];
   double beta = par[2];
   auto exp = TMath::Exp(-(t-t_shift)/alpha);
   double res = (1.0-exp)*(1.0+beta*t);
   //std::cout << res << "  res\n";
   //double integral = integral_acc(Params::tMax, par)-integral_acc(Params::tMin, par);
   if (res <= 0.0) return 0.0;
   return res;
   //return (1.0+beta*t);   
}
double pdfs::integral_acc(double t, const double *par){
   double t_shift = par[0];
   double alpha = par[1];
   double beta = par[2];
   auto exp = TMath::Exp(-(t-t_shift)/alpha);
     double res = (1.0-exp);
     res*= 1.0+beta*t;
   if (res <= 0.0) return 0.0;

   //return 1./2.*beta*t*t+t +alpha*TMath::Exp((t_shift-t)/alpha)*(1.0+beta+alpha*beta);
   return TMath::Power((t - t_shift),alpha)/(1.0+beta*TMath::Power((t -t_shift), alpha));
}
double pdfs::func1D_f(double *x, const double *par) // 1D untagged function of decay to fermion final state
{
  auto eta = par[3];
  //return Params::N*((1.0+eta)*pdfs::B_to_f(x, par) + (1.0-eta)*pdfs::Bbar_to_f(x, par));
  //auto t_min = Params::tMin;	
  //auto t_max = Params::tMax;	
  //auto integral = integral_func1D(&t_max, par, 1.)- integral_func1D(&t_min, par, 1.);
  return (pdfs::B_to_f(x, par) + eta*pdfs::Bbar_to_f(x, par));
}
double pdfs::func1D_fbar(double *x, const double *par)//1D untagged function of decay to antifermion final state
{
  auto eta = par[3];
  //return Params::Nbar*((1.0+eta)*pdfs::B_to_fbar(x, par) + (1.0-eta)*pdfs::Bbar_to_fbar(x, par));
  //auto t_min = Params::tMin;	
  //auto t_max = Params::tMax;	
  //auto integral = integral_func1D(&t_max, par, -1.)- integral_func1D(&t_min, par, -1.);
  return (1.0/eta*pdfs::B_to_fbar(x, par) + pdfs::Bbar_to_fbar(x, par));
}
double pdfs::integral_func1D(double *x, const double *par, double m){
  auto t = x[0];
  auto R = par[0];//*(-2.0)*Params::dM/Params::dGamma;
  auto I = par[1];
  auto Gamma = par[2];
  auto eta = par[3];
  if (m<0) eta = 1./eta;
  auto dG = Params::dGamma;	
  auto dM = Params::dM;

  auto res = -TMath::Exp(-Gamma*t-dG*t/2)*(sqrt(R*R*R*R+(2.*I*I-2.)*R*R+I*I*I*I+2.*I*I+1.)
	*((dG*dG-4.*Gamma*Gamma)*dM*eta*TMath::Exp(dG*t/2.)*TMath::Sin(dM*t)
	+(4.*Gamma*Gamma*Gamma-Gamma*dG*dG)*eta*TMath::Exp(dG*t/2.)*TMath::Cos(dM*t)
	+(((-dG)-2.*Gamma)*dM*dM-Gamma*Gamma*dG-2.*Gamma*Gamma*Gamma)*eta*TMath::Exp(dG*t)
	+((dG-2.*Gamma)*dM*dM+Gamma*Gamma*dG-2.*Gamma*Gamma*Gamma)*eta)
	+(((R*R+I*I-1.)*dG*dG-4.*Gamma*Gamma*R*R-4.*Gamma*Gamma*I*I+4.*Gamma*Gamma)*dM
	-m*2.*Gamma*I*dG*dG+m*8.*Gamma*Gamma*Gamma*I)*TMath::Exp(dG*t/2.)*TMath::Sin(dM*t)
	+(m*(8.*Gamma*Gamma*I-2.*I*dG*dG)*dM+((-Gamma*R*R)-Gamma*I*I+Gamma)*dG*dG
	+4.*Gamma*Gamma*Gamma*R*R+4*Gamma*Gamma*Gamma*I*I-4.*Gamma*Gamma*Gamma)*TMath::Exp(dG*t/2.)*TMath::Cos(dM*t)
	+((((-R*R)+m*2.*R-I*I-1.)*dG-2.*Gamma*R*R+m*4.*Gamma*R-2.*Gamma*I*I
	-2.*Gamma)*dM*dM+((-Gamma*Gamma*R*R)+m*2.*Gamma*Gamma*R-Gamma*Gamma*I*I-Gamma*Gamma)*dG
	-2.*Gamma*Gamma*Gamma*R*R+m*4.*Gamma*Gamma*Gamma*R -2.*Gamma*Gamma *Gamma *I*I-2.*Gamma*Gamma*Gamma)*TMath::Exp(dG*t)
	+((R*R+2.*m*R+I*I+1.)*dG-2.*Gamma*R*R-m*4.*Gamma*R-2.*Gamma*I*I-2.*Gamma)
	*dM*dM+(Gamma*Gamma*R*R+m*2.*Gamma*Gamma*R+Gamma*Gamma*I*I+Gamma*Gamma)*dG-2.*Gamma*Gamma*Gamma*R*R
	-m*4.*Gamma*Gamma*Gamma*R-2.*Gamma*Gamma*Gamma*I*I-2.*Gamma*Gamma*Gamma);

  res/=((dG*dG-4*Gamma*Gamma)*dM*dM+Gamma*Gamma*dG*dG-4.0*Gamma*Gamma*Gamma*Gamma);
  return res;
}
  	  
double pdfs::func1D_f_simple(double *x, double *par){
  auto t = x[0];
  auto imz = par[1];
  auto Gamma = par[2];
  auto eta = par[3];
  auto exp = TMath:: Exp(-Gamma*t);
  return exp*(1.+eta+(1.-eta)*TMath::Cos(Params::dM*t)-2.*imz*TMath::Sin(Params::dM*t));
}
double pdfs::func1D_fbar_simple(double *x, double *par){
  auto t = x[0];
  auto imz = par[1];
  auto Gamma = par[2];
  auto eta = par[3];
  auto exp = TMath:: Exp(-Gamma*t);
  return exp*(1.+1./eta+(1.-1./eta)*TMath::Cos(Params::dM*t)+2.*imz*TMath::Sin(Params::dM*t));
}
//////////////// Normalisations /////////////////////////////////////////////
/*
double pdfs::normalization_B_to_f(double rez, double imz, double m)
{
  double dGamma = Params::dGamma;
  double Gamma = Params::Gamma;
  double dM = Params::dM;
  double t_min = global.tMin_gen;
  double t_max = global.tMax_gen;

  double correction = dGamma*dGamma + Gamma*Gamma;

  double e_hypcosine_top = 2.*dGamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::SinH(dGamma*t_max*0.5)/correction;
  e_hypcosine_top += -2.*Gamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::CosH(dGamma*t_max*0.5)/correction; 
  double e_hypcosine_bottom = 2.*dGamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::SinH(dGamma*t_min*0.5)/correction;
  e_hypcosine_bottom += -2.*Gamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::CosH(dGamma*t_min*0.5)/correction;

  correction = Gamma*Gamma + 4.*dM*dM;
  double e_cosine_top = 4.*dM*TMath::Exp(-Gamma*t_max*0.5)*TMath::Sin(dM*t_max)/correction;
  e_cosine_top += -2.*Gamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::Cos(dM*t_max)/correction;
  double e_cosine_bottom = 4.*dM*TMath::Exp(-Gamma*t_min*0.5)*TMath::Sin(dM*t_min)/correction;
  e_cosine_bottom += -2.*Gamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::Cos(dM*t_min)/correction;

  correction = dGamma*dGamma + Gamma*Gamma;
  double e_hypsine_top = -2.*Gamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::SinH(dGamma*t_max*0.5)/correction;
  e_hypsine_top += -2.*dGamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::CosH(dGamma*t_max*0.5)/correction;
  double e_hypsine_bottom = -2.*Gamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::SinH(dGamma*t_min*0.5)/correction;
  e_hypsine_bottom += -2.*dGamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::CosH(dGamma*t_min*0.5)/correction;

  correction = Gamma*Gamma + 4.*dM*dM;
  double e_sine_top = -2.*Gamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::Sin(dM*t_max)/correction;
  e_sine_top += -4.*dM*TMath::Exp(-Gamma*t_max*0.5)*TMath::Cos(dM*t_max)/correction;
  double e_sine_bottom = -2.*Gamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::Sin(dM*t_min)/correction;
  e_sine_bottom += -4.*dM*TMath::Exp(-Gamma*t_min*0.5)*TMath::Cos(dM*t_min)/correction;

  double z_sq = rez*rez+imz*imz;
  double N0 = (1+z_sq)*(e_hypcosine_top-e_hypcosine_bottom)+(1-z_sq)*(e_cosine_top-e_cosine_bottom);
  N0+= m*(-2.*rez*(e_hypsine_top-e_hypsine_bottom)-2.*imz*(e_sine_top-e_sine_bottom));
  return N0;
  return N0*global.nbins/(t_max-t_min);
}
double pdfs::normalization_CP(double rez, double imz)
{
  double dGamma = Params::dGamma;
  double Gamma = Params::Gamma;
  double dM = Params::dM;
  double t_min = global.tMin_gen;
  double t_max = global.tMax_gen;

  double correction = dGamma*dGamma + Gamma*Gamma;
  double e_hypcosine_top = 2.*dGamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::SinH(dGamma*t_max*0.5)/correction;
  e_hypcosine_top += -2.*Gamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::CosH(dGamma*t_max*0.5)/correction;
  double e_hypcosine_bottom = 2.*dGamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::SinH(dGamma*t_min*0.5)/correction;
  e_hypcosine_bottom += -2.*Gamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::CosH(dGamma*t_min*0.5)/correction;


  correction = Gamma*Gamma + 4.*dM*dM;
  double e_cosine_top = 4.*dM*TMath::Exp(-Gamma*t_max*0.5)*TMath::Sin(dM*t_max)/correction;
  e_cosine_top += -2.*Gamma*TMath::Exp(-Gamma*t_max*0.5)*TMath::Cos(dM*t_max)/correction;
  double e_cosine_bottom = 4.*dM*TMath::Exp(-Gamma*t_min*0.5)*TMath::Sin(dM*t_min)/correction;
  e_cosine_bottom += -2.*Gamma*TMath::Exp(-Gamma*t_min*0.5)*TMath::Cos(dM*t_min)/correction;

  auto zeta_CP = TMath::Sqrt((1-rez*rez+imz*imz)*(1-rez*rez+imz*imz)+4.*rez*rez*imz*imz);
  auto CP = Params::omega2;
  double N0 = zeta_CP*(e_hypcosine_top-e_hypcosine_bottom-e_cosine_top+e_cosine_bottom);
  return N0;
  return N0*global.nbins/(t_max-t_min);
}
*/

}// cpt_b0_analysis namespace
