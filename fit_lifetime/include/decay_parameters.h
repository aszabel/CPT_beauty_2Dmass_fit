
#ifndef DECAY_PARAMETERS_H
#define DECAY_PARAMETERS_H


#include <stdio.h>


namespace cpt_b0_analysis
{
	/*
    class Global{
  
    static double tMin = 1.0;
    static double tMax = 15.0;
    static double tMin_gen = 0.0;
    static double tMax_gen = 18.0;
    static double Mmin = 1825.; 
    static double Mmax = 1920.;
    
    static int N = 1000000;
    static int Nbar = 1000000; 
    static int nbins = 10;
    static int nbins1D = 200;
    static int npar = 37;

    //double bkg_sb_percent =0.064;
    //double bkg_sb_percent =0.1670;
    static double bkg_sb_percent = 0.18;//0.177875;
    static double bkg_sb_percent_bar =0.18;//0.177711;
    static double signal_percent = 1.0-0.162;
    static double signal_percent_bar = signal_percent;
    //double bkg_sb_percent =1.0;
    //double signal_percent =0.0;
    static double bkg_Bu_percent = (1.0-bkg_sb_percent)*(1.0-signal_percent);
    

    static bool eq_bins = false;
    static bool simfit = true;
    static bool sb_bkg = bkg_sb_percent != 0.0; 
    static bool bu_bkg = bkg_Bu_percent != 0.0;
    static bool b0_sig = signal_percent*(1.0-bkg_sb_percent) != 0.0;

    static double seed = 0; /// seed used by the MC generator
    static bool data_from_file = true;
    static bool r_hist_from_file = true;
    static bool chi2 = false;
  };*/
  class Params {
	  public:
    Params() {}
    /*
    Params(constexpr double Gamma_, constexpr double dGamma_, constexpr double dM_, double Rez_, double Imz_, double omega2_,
           double eta_, double A_D_, double A_dir_, 
           double sigma_, 
           int N_, int Nbar_, int Nr_, int Nrbins_, 
           int npar_,
           double alpha_, double beta_, double t_shift_,
           bool time_acc_, bool conv_):
           Gamma(Gamma_), dGamma(dGamma_), dM(dM_), Rez(Rez_), Imz(Imz_), omega2(omega2_),
           eta(eta_), A_D(A_D_), A_dir(A_dir_),
           sigma(sigma_), 
           N(N_), Nbar(Nbar_), Nr(Nr_), Nrbins(Nrbins_),
           npar(npar_),
           alpha(alpha_), beta(beta_), t_shift(t_shift_),
           time_acc(time_acc_), conv(conv_)
           {}
*/

    static constexpr double tMin = 1.0;
    static constexpr double tMax = 15.0;
    static constexpr double tMin_gen = 0.0;
    static constexpr double tMax_gen = 18.0;

    static constexpr double Gamma = 658.3/1000.0;
    static constexpr double dGamma =-0.002*Gamma;
    static constexpr double dM =0.77*Gamma;
    static constexpr double Rez = 0.0;
    static constexpr double Imz = 0.0;//0.01;
    static constexpr double omega2 = 1.0;//0.9998;
    static constexpr double eta = 1.00246187;
    static constexpr double A_D = 0.0;
    static constexpr double A_dir = 0.0;
    static constexpr double sigma = 0.00045;
    static constexpr int N = 0;//int(((1.0 - global.bkg_sb_percent)*global.signal_percent)*double(global.N));
    static constexpr int Nr= 10000;
    static constexpr int Nbar = 0;//int(((1.0-global.bkg_sb_percent)*global.signal_percent)*double(global.Nbar));
    static constexpr int Nrbins = 50;
    static constexpr int npar = 10;
    static constexpr bool conv = true;

    static constexpr bool time_acc = true;
    //double alpha = 0.000554;
    //double beta = -36.;
    //double t_shift = 0.000197;
    //double alpha = 0.0004;
    //double beta = -40.;
    //double t_shift = 0.0009;
    static constexpr double alpha = 0.24;
    static constexpr double beta = -0.0385;
    static constexpr double t_shift = 0.03;//0.00006;

//    bool _first = true;
  };
/*
  struct DmassParams{
    constexpr int npar = 12;
    double sigma = 10.8939;
    double mean = 1870.04;
    double exp_slope = 0.00500569;
    double s12 = 0.58429*sigma;
    double f12 = 0.381135	;
    double fcomb = 0.18;//0.177875;
    double fcomb_bar = 0.18;//0.177711;
    double fcomb_gen = fcomb_bar;
    double fcomb_bar_gen = fcomb;
  };
     
*/


/*  Params b0_pars; 

  Params combBKG(
  2.0*673.0, 0.0, 0.510, 0.0, 0.0, 1.0,
  -0.008, 0.0, 0.0, 
  0.0, 
  int(global.bkg_sb_percent*global.N), int(global.bkg_sb_percent*global.Nbar), 0, 1,
  3,
  0.0, 0.0, 0.0,
  false, false);

  Params BuBKG(
  2.0*610.5, 0.0, 0.0, 0.0, 0.0, 1.0,
  -0.0006, 0.0, 0.0,
  b0_pars.sigma,
  int(global.bkg_Bu_percent*b0_pars.N), int(global.bkg_Bu_percent*b0_pars.Nbar), b0_pars.Nr, b0_pars.Nrbins,
  3,                                                                               //acceptance the same as for singnal 
  b0_pars.alpha, b0_pars.beta, b0_pars.t_shift,
  true, true);
*/    
  

}
#endif
