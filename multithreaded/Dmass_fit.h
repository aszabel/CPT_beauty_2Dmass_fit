/**
  * Adam Szabelski 28.01.2020
  *
  * @Dmass_fit.h - declaration of the Dmass fitting class
  *
**/
#ifndef DMASS_FIT_H
#define DMASS_FIT_H

#include "TH2D.h"
#include<stdio.h>
#include "M_B_2missPT_fit.h"

namespace cpt_b0_analysis{
const int nprm=9;
const int ncontr=6;
const int nsigBu = 6;
const int nBs = 4;
const int nB0D = 7;
const int nBuD = 7;
const int nvar = 87;
class Dmass{
  public:
    Dmass(std::vector<std::pair<double, double>>& vect_Dmass, TH2D **hcorr, double minx, double maxx, double miny, double maxy, bool sign = true);
    ~Dmass();
    void SetExternal(double xx[ncontr][nprm], double dxx[ncontr][nprm], double yy[nsigBu+nBs+nB0D+nBuD], double dyy[nsigBu+nBs+nB0D+nBuD]);
    double operator()(const double *par);

    double Dmass_func(double *x, const double *par);
    double gausPDF(double *x, const double *par, int bin_x, int bin_y);
    double gausPDF_B0D(double *x, const double *par, int bin_x, int bin_y);
    double gausPDF_BuD(double *x, const double *par, int bin_x, int bin_y);
    double crystalballPDF(double *x, const double *par, int bin_x, int bin_y);
    double crystalballPDF2(double *x, const double *par, int bin_x, int bin_y);
    double DoubleSidedCrystalballFunction_B0D(double *x, const double *par, int bin_x, int bin_y);
    double DoubleSidedCrystalballFunction_BuD(double *x, const double *par, int bin_x, int bin_y);
    double exponentPDF(double *x, const double *par, int bin_x, int bin_y);
    double chebyshevPDF(double *x, const double *par, int bin_x, int bin_y);
    double calcChi2(const double *param);
    void Draw(const double *param, const char *name, double CovMTx[nvar*nvar]);
    double gaussian_error(double mean, double error, const double *param, int ipar);
    void generate (const double *param, int N);
    void normalize_histogram(TH2D &hist);
    void GetCorrBin(int &bin_x, int &bin_y, double m_rec, double m_corr);

    std::vector<std::pair<double, double>> vect_Dmass;
    TH2D hdata;
    TH2D **hcorr;
    bool sign;
    int nbins;
    double nevents; 
    double bin_width, bin_widthy;
    double tabcorr[10][10][100];
    double minx, miny, maxx, maxy, xx[ncontr][nprm], dxx[ncontr][nprm], yy[nsigBu+nBs+nB0D+nBuD], dyy[nsigBu+nBs+nB0D+nBuD];
    //kde_fit kde_graphs;
    mb_2misspt_fit *mb_2misspt_graphs;

};
class Dmass_autofit: public Dmass{
  public:
    Dmass_autofit (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_autofit(){};
    double operator()(double *x, double *par);   
};
class Dmass_genExpo: public Dmass{
  public:
    Dmass_genExpo (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genExpo(){};
    double operator()(double *x, double *par);   
};
class Dmass_genCheb: public Dmass{
  public:
    Dmass_genCheb (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genCheb(){};
    double operator()(double *x, double *par);   
};
class Dmass_genGaus: public Dmass{
  public:
    Dmass_genGaus (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genGaus(){};
    double operator()(double *x, double *par);   
};
class Dmass_genGausD: public Dmass{
  public:
    Dmass_genGausD (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genGausD(){};
    double operator()(double *x, double *par);   
};
class Dmass_genD0B: public Dmass{
  public:
    Dmass_genD0B (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genD0B(){};
    double operator()(double *x, double *par);
};

class Dmass_genDuB: public Dmass{
  public:
    Dmass_genDuB (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genDuB(){};
    double operator()(double *x, double *par);
};

class Dmass_genGaus2: public Dmass{
  public:
    Dmass_genGaus2 (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genGaus2(){};
    double operator()(double *x, double *par);   
};

class Dmass_genCB2: public Dmass{
  public:
    Dmass_genCB2 (std::vector<std::pair<double, double>> vect_Dmass_, TH2D **hcorr_, double minx_, double maxx_, double miny_, double maxy_);
    ~Dmass_genCB2(){};
    double operator()(double *x, double *par);   
};

    
}
#endif 
