/**
 *  Adam Szabelski 13.06.2019 tool for generation and analysis of CPT conservation in semileptonic decays 
 *  A list of semileptonic decay time functions to be fitted and or generated
 *
 *
 *  @file pdfs_cpt_model.h
*/
#ifndef PDFS_CPT_MODEL_H
#define PDFS_CPT_MODEL_H

#include "decay_parameters.h"

namespace cpt_b0_analysis{

class pdfs{

  public:
    
    pdfs(){}


    double B_to_f(double *x, const double *par);

    double Bbar_to_f(double *x, const double *par);

    double calc_B_to_f(double *x, const double *par, double m);

    double B_to_fbar(double *x, const double *par);

    double Bbar_to_fbar(double *x, const double *par);

    double calcB_to_fbar(double *x, const double *par, double m);

    double fitf_unasym(double *x, double *par);

    double fitf_decay(double *x, double *par);

    double fitf_simult(double *x, double *par);

    double time_acceptance(double *x, const double *par);

    double func1D_f(double *x, const double *par);

    double func1D_fbar(double *x, const double *par);

    double func1D_f_simple(double *x, double *par);

    double func1D_fbar_simple(double *x, double *par);

//    double normalization_B_to_f(double rez, double imz, double m);

//    double normalization_CP(double rez, double imz);
    double integral_func1D(double *x, const double *par, double m);
    double integral_acc(double y, const double *par);
};
}// cpt_b0_analysis namspace

#endif /*PDFS_CPT_MODEL_H */
