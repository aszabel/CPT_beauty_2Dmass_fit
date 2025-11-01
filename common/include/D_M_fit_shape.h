#ifndef DM_FIT_H
#define DM_FIT_H
#include <vector>

#include "BasicShapes.h"
#include "ChebyshevPDF.h"
#include "pdf_interface.h"

namespace cpt_b0_analysis {

class DoubleSidedCrystalballPlusGaussPDF : public PDFInterface {
public:
	DoubleSidedCrystalballPlusGaussPDF() {};
	DoubleSidedCrystalballPlusGaussPDF(const DoubleSidedCrystalballPlusGaussPDF &other) = default;
	DoubleSidedCrystalballPlusGaussPDF(DoubleSidedCrystalballPlusGaussPDF &&) = default;
	~DoubleSidedCrystalballPlusGaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	GaussPDF Gauss;
	DoubleSidedCrystalBallPDF DSCB;
};

class DoubleSidedCrystalballPlusExpPDF : public PDFInterface {
public:
	DoubleSidedCrystalballPlusExpPDF() {};
	DoubleSidedCrystalballPlusExpPDF(const DoubleSidedCrystalballPlusExpPDF &other) = default;
	DoubleSidedCrystalballPlusExpPDF(DoubleSidedCrystalballPlusExpPDF &&) = default;
	~DoubleSidedCrystalballPlusExpPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	ExponentPDF Exp;
	DoubleSidedCrystalBallPDF DSCB;
};

class JohnsonPlusGaussPDF : public PDFInterface {
public:
	JohnsonPlusGaussPDF() {};
	JohnsonPlusGaussPDF(const JohnsonPlusGaussPDF &other) = default;
	JohnsonPlusGaussPDF(JohnsonPlusGaussPDF &&) = default;
	~JohnsonPlusGaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	GaussPDF Gauss;
	JohnsonPDF JSU;
};

class JohnsonPlusExpPDF : public PDFInterface {
public:
	JohnsonPlusExpPDF() {};
	JohnsonPlusExpPDF(const JohnsonPlusExpPDF &other) = default;
	JohnsonPlusExpPDF(JohnsonPlusExpPDF &&) = default;
	~JohnsonPlusExpPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	ExponentPDF Exp;
	JohnsonPDF JSU;
};

class JohnsonPlusDoubleSidedCrystalBallPDF : public PDFInterface {
public:
	JohnsonPlusDoubleSidedCrystalBallPDF() {};
	JohnsonPlusDoubleSidedCrystalBallPDF(const JohnsonPlusDoubleSidedCrystalBallPDF &other) =
		default;
	JohnsonPlusDoubleSidedCrystalBallPDF(JohnsonPlusDoubleSidedCrystalBallPDF &&) = default;
	~JohnsonPlusDoubleSidedCrystalBallPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	DoubleSidedCrystalBallPDF DSCB;
	JohnsonPDF JSU;
};
}  // namespace cpt_b0_analysis

#endif
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
