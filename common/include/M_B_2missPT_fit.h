#ifndef M_B_PDF_H
#define M_B_PDF_H
#include <functional>
#include <vector>

#include "BasicShapes.h"
#include "pdf_interface.h"

namespace cpt_b0_analysis {

class RaisedCosinePlusGaussPDF : public PDFInterface {
public:
	RaisedCosinePlusGaussPDF() {};
	RaisedCosinePlusGaussPDF(const RaisedCosinePlusGaussPDF &other) = default;
	RaisedCosinePlusGaussPDF(RaisedCosinePlusGaussPDF &&) = default;
	~RaisedCosinePlusGaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 3; };	// RaisedCosine + 2 x Gauss
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	RaisedCosinePDF Cos;
	DoubleGaussPDF DGauss;
};

class SkewNormalPlusCBPDF : public PDFInterface {
public:
	SkewNormalPlusCBPDF() {};
	SkewNormalPlusCBPDF(const SkewNormalPlusCBPDF &other) = default;
	SkewNormalPlusCBPDF(SkewNormalPlusCBPDF &&) = default;
	~SkewNormalPlusCBPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };	// SkewNormal + CB
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	SkewNormalPDF SkewNorm;
	CrystalBallPDF CB;
};

class JohnsonPlusCBPDF : public PDFInterface {
public:
	JohnsonPlusCBPDF() {};
	JohnsonPlusCBPDF(const JohnsonPlusCBPDF &other) = default;
	JohnsonPlusCBPDF(JohnsonPlusCBPDF &&) = default;
	~JohnsonPlusCBPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };	// Johnson + CB
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	JohnsonPDF JSU;
	CrystalBallPDF CB;
};

class SkewNormalPlusCBPlusDoubleGaussPDF : public PDFInterface {
public:
	SkewNormalPlusCBPlusDoubleGaussPDF() {};
	SkewNormalPlusCBPlusDoubleGaussPDF(const SkewNormalPlusCBPlusDoubleGaussPDF &other) = default;
	SkewNormalPlusCBPlusDoubleGaussPDF(SkewNormalPlusCBPlusDoubleGaussPDF &&) = default;
	~SkewNormalPlusCBPlusDoubleGaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 4; };	// SkewNormal + CB
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	SkewNormalPDF SkewNorm;
	CrystalBallPDF CB;
	DoubleGaussPDF Gauss;
};

class JohnsonPlusCBPlusDoubleGaussPDF : public PDFInterface {
public:
	JohnsonPlusCBPlusDoubleGaussPDF() {};
	JohnsonPlusCBPlusDoubleGaussPDF(const JohnsonPlusCBPlusDoubleGaussPDF &other) = default;
	JohnsonPlusCBPlusDoubleGaussPDF(JohnsonPlusCBPlusDoubleGaussPDF &&) = default;
	~JohnsonPlusCBPlusDoubleGaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 4; };	// Johnson + CB
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	JohnsonPDF JSU;
	CrystalBallPDF CB;
	DoubleGaussPDF Gauss;
};
}  // namespace cpt_b0_analysis

#endif
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
