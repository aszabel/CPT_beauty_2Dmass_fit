#include "M_B_2missPT_fit.h"

#include <Math/DistFunc.h>	// normal_cdf

#include <boost/math/special_functions/owens_t.hpp>

#include "Math/GaussIntegrator.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/WrappedTF1.h"
#include "TF1.h"
#include "TMath.h"

namespace cpt_b0_analysis {

double RaisedCosinePlusGaussPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double f = abs(par[0]);
	double par1[2];
	double par2[5];
	par1[0] = par[1];  // mean_cos
	par1[1] = par[2];  // sigma_cos
	par2[0] = par[3];  // f_gauss
	par2[1] = par[4];  // mean_1
	par2[2] = par[5];  // sigma_1
	par2[3] = par[6];  // mean_2
	par2[4] = par[7];  // sigma_2

	if (component == 0)
		return (1 - f) * Cos.EvalPDF(xx, par1);
	else if (component == 1)
		return f * DGauss.EvalPDF(xx, par2, 0);
	else if (component == 2)
		return f * DGauss.EvalPDF(xx, par2, 1);
	else
		return (1 - f) * Cos.EvalPDF(xx, par1) + f * DGauss.EvalPDF(xx, par2);
}

void RaisedCosinePlusGaussPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[2];
	double par2[5];
	par1[0] = par[1];  // mean_cos
	par1[1] = par[2];  // sigma_cos
	par2[0] = par[3];  // f_gauss
	par2[1] = par[4];  // mean_1
	par2[2] = par[5];  // sigma_1
	par2[3] = par[6];  // mean_2
	par2[4] = par[7];  // sigma_2

	Cos.CalcIntegral(par1, min, max);
	DGauss.CalcIntegral(par2, min, max);
}

double SkewNormalPlusCBPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double f = abs(par[0]);
	double par1[3];
	double par2[4];
	par1[0] = par[1];  // mean
	par1[1] = par[2];  // sigma
	par1[2] = par[3];  // skew
	par2[0] = par[4];  // meanCB
	par2[1] = par[5];  // sigmaCB
	par2[2] = par[6];  // alpha
	par2[3] = par[7];  // n

	if (component == 0)
		return (1 - f) * SkewNorm.EvalPDF(xx, par1);
	else if (component == 1)
		return f * CB.EvalPDF(xx, par2);
	else
		return (1 - f) * SkewNorm.EvalPDF(xx, par1) + f * CB.EvalPDF(xx, par2);
}

void SkewNormalPlusCBPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[3];
	double par2[4];
	par1[0] = par[1];  // mean
	par1[1] = par[2];  // sigma
	par1[2] = par[3];  // skew
	par2[0] = par[4];  // meanCB
	par2[1] = par[5];  // sigmaCB
	par2[2] = par[6];  // alpha
	par2[3] = par[7];  // n

	SkewNorm.CalcIntegral(par1, min, max);
	CB.CalcIntegral(par2, min, max);
}

double JohnsonPlusCBPDF::EvalPDF(const double *xx, const double *par, const int component) {
	double f = abs(par[0]);
	double par1[4];
	double par2[4];
	par1[0] = par[1];  // xi
	par1[1] = par[2];  // lambda
	par1[2] = par[3];  // gamma
	par1[3] = par[4];  // delta
	par2[0] = par[5];  // mean
	par2[1] = par[6];  // sigma
	par2[2] = par[7];  // alpha
	par2[3] = par[8];  // n

	if (component == 0)
		return (1 - f) * JSU.EvalPDF(xx, par1);
	else if (component == 1)
		return f * CB.EvalPDF(xx, par2);
	else
		return (1 - f) * JSU.EvalPDF(xx, par1) + f * CB.EvalPDF(xx, par2);
}

void JohnsonPlusCBPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[4];
	double par2[4];
	par1[0] = par[1];  // xi
	par1[1] = par[2];  // lambda
	par1[2] = par[3];  // gamma
	par1[3] = par[4];  // delta
	par2[0] = par[5];  // mean
	par2[1] = par[6];  // sigma
	par2[2] = par[7];  // alpha
	par2[3] = par[8];  // n

	JSU.CalcIntegral(par1, min, max);
	CB.CalcIntegral(par2, min, max);
}

double SkewNormalPlusCBPlusDoubleGaussPDF::EvalPDF(const double *xx, const double *par,
												   const int component) {
	double f_gauss = abs(par[0]);
	double f_12 = abs(par[1]);
	double par1[3];
	double par2[4];
	double par3[5];
	par1[0] = par[2];	// mean
	par1[1] = par[3];	// sigma
	par1[2] = par[4];	// skew
	par2[0] = par[5];	// meanCB
	par2[1] = par[6];	// sigmaCB
	par2[2] = par[7];	// alpha
	par2[3] = par[8];	// n
	par3[0] = par[9];	// f gauss12
	par3[1] = par[10];	// mean gauss 1
	par3[2] = par[11];	// sigma gauss 1
	par3[3] = par[12];	// mean gauss 2
	par3[4] = par[13];	// sigma gauss 2

	if (component == 0)
		return (1.0 - f_gauss) * (1.0 - f_12) * SkewNorm.EvalPDF(xx, par1);
	else if (component == 1)
		return (1.0 - f_gauss) * f_12 * CB.EvalPDF(xx, par2);
	else if (component == 2)
		return f_gauss * Gauss.EvalPDF(xx, par3, 0);
	else if (component == 3)
		return f_gauss * Gauss.EvalPDF(xx, par3, 1);
	else
		return (1 - f_gauss) * (1.0 - f_12) * SkewNorm.EvalPDF(xx, par1) +
			   (1 - f_gauss) * f_12 * CB.EvalPDF(xx, par2) + f_gauss * Gauss.EvalPDF(xx, par3);
}

void SkewNormalPlusCBPlusDoubleGaussPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[3];
	double par2[4];
	double par3[5];
	par1[0] = par[2];	// mean
	par1[1] = par[3];	// sigma
	par1[2] = par[4];	// skew
	par2[0] = par[5];	// meanCB
	par2[1] = par[6];	// sigmaCB
	par2[2] = par[7];	// alpha
	par2[3] = par[8];	// n
	par3[0] = par[9];	// f gauss12
	par3[1] = par[10];	// mean gauss 1
	par3[2] = par[11];	// sigma gauss 1
	par3[3] = par[12];	// mean gauss 2
	par3[4] = par[13];	// sigma gauss 2

	SkewNorm.CalcIntegral(par1, min, max);
	CB.CalcIntegral(par2, min, max);
	Gauss.CalcIntegral(par3, min, max);
}

double JohnsonPlusCBPlusDoubleGaussPDF::EvalPDF(const double *xx, const double *par,
												const int component) {
	double f_gauss = abs(par[0]);
	double f_12 = abs(par[1]);
	double par1[4];
	double par2[4];
	double par3[5];
	par1[0] = par[2];	// xi
	par1[1] = par[3];	// lambda
	par1[2] = par[4];	// gamma
	par1[3] = par[5];	// delta
	par2[0] = par[6];	// mean
	par2[1] = par[7];	// sigma
	par2[2] = par[8];	// alpha
	par2[3] = par[9];	// n
	par3[0] = par[10];	// f gauss12
	par3[1] = par[11];	// mean gauss 1
	par3[2] = par[12];	// sigma gauss 1
	par3[3] = par[13];	// mean gauss 2
	par3[4] = par[14];	// sigma gauss 2

	if (component == 0)
		return (1.0 - f_gauss) * (1.0 - f_12) * JSU.EvalPDF(xx, par1);
	else if (component == 1)
		return (1.0 - f_gauss) * f_12 * CB.EvalPDF(xx, par2);
	else if (component == 2)
		return f_gauss * Gauss.EvalPDF(xx, par3, 0);
	else if (component == 3)
		return f_gauss * Gauss.EvalPDF(xx, par3, 1);
	else
		return (1.0 - f_gauss) * (1.0 - f_12) * JSU.EvalPDF(xx, par1) +
			   (1.0 - f_gauss) * f_12 * CB.EvalPDF(xx, par2) + f_gauss * Gauss.EvalPDF(xx, par3);
}

void JohnsonPlusCBPlusDoubleGaussPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[4];
	double par2[4];
	double par3[5];
	par1[0] = par[2];	// xi
	par1[1] = par[3];	// lambda
	par1[2] = par[4];	// gamma
	par1[3] = par[5];	// delta
	par2[0] = par[6];	// mean
	par2[1] = par[7];	// sigma
	par2[2] = par[8];	// alpha
	par2[3] = par[9];	// n
	par3[0] = par[10];	// f gauss12
	par3[1] = par[11];	// mean gauss 1
	par3[2] = par[12];	// sigma gauss 1
	par3[3] = par[13];	// mean gauss 2
	par3[4] = par[14];	// sigma gauss 2

	JSU.CalcIntegral(par1, min, max);
	CB.CalcIntegral(par2, min, max);
	Gauss.CalcIntegral(par3, min, max);
}

double SidebandsPDF::EvalPDF(const double *xx, const double *par,
												const int component) {
	double f_12 = abs(par[0]);
	double par1[15];
	double par2[15];
	par1[0] = par[1];	// f_gauss
	par1[1] = par[2];	// f12
	par1[2] = par[3];	// xi
	par1[3] = par[4];	// lambda
	par1[4] = par[5];	// gamma
	par1[5] = par[6];	// delta
	par1[6] = par[7];	// mean
	par1[7] = par[8];	// sigma
	par1[8] = par[9];	// alpha
	par1[9] = par[10];	// n
	par1[10] = par[11];	// f gauss12
	par1[11] = par[12];	// mean gauss 1
	par1[12] = par[13];	// sigma gauss 1
	par1[13] = par[14];	// mean gauss 2
	par1[14] = par[15];	// sigma gauss 2

	par2[0] = par[16];	// f_gauss
	par2[1] = par[17];	// f12
	par2[2] = par[18];	// xi
	par2[3] = par[19];	// lambda
	par2[4] = par[20];	// gamma
	par2[5] = par[21];	// delta
	par2[6] = par[22];	// mean
	par2[7] = par[23];	// sigma
	par2[8] = par[24];	// alpha
	par2[9] = par[25];	// n
	par2[10] = par[26];	// f gauss12
	par2[11] = par[27];	// mean gauss 1
	par2[12] = par[28];	// sigma gauss 1
	par2[13] = par[29];	// mean gauss 2
	par2[14] = par[30];	// sigma gauss 2

	if (component == 0)
		return (1.0 - f_12) * Left.EvalPDF(xx, par1);
	else if (component == 1)
		return f_12 * Right.EvalPDF(xx, par2);
	else
		return (1.0 - f_12) * Left.EvalPDF(xx, par1) +
			   f_12 * Right.EvalPDF(xx, par2);
}

void SidebandsPDF::CalcIntegral(const double *par, double min, double max) {
	double par1[15];
	double par2[15];
	par1[0] = par[1];	// f_gauss
	par1[1] = par[2];	// f12
	par1[2] = par[3];	// xi
	par1[3] = par[4];	// lambda
	par1[4] = par[5];	// gamma
	par1[5] = par[6];	// delta
	par1[6] = par[7];	// mean
	par1[7] = par[8];	// sigma
	par1[8] = par[9];	// alpha
	par1[9] = par[10];	// n
	par1[10] = par[11];	// f gauss12
	par1[11] = par[12];	// mean gauss 1
	par1[12] = par[13];	// sigma gauss 1
	par1[13] = par[14];	// mean gauss 2
	par1[14] = par[15];	// sigma gauss 2

	par2[0] = par[16];	// f_gauss
	par2[1] = par[17];	// f12
	par2[2] = par[18];	// xi
	par2[3] = par[19];	// lambda
	par2[4] = par[20];	// gamma
	par2[5] = par[21];	// delta
	par2[6] = par[22];	// mean
	par2[7] = par[23];	// sigma
	par2[8] = par[24];	// alpha
	par2[9] = par[25];	// n
	par2[10] = par[26];	// f gauss12
	par2[11] = par[27];	// mean gauss 1
	par2[12] = par[28];	// sigma gauss 1
	par2[13] = par[29];	// mean gauss 2
	par2[14] = par[30];	// sigma gauss 2

	Left.CalcIntegral(par1, min, max);
	Right.CalcIntegral(par2, min, max);
}
}  // namespace cpt_b0_analysis

// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
