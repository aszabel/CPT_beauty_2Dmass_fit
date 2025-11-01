#ifndef BASIC_SHAPES_H
#define BASIC_SHAPES_H
#include <functional>
#include <vector>

#include "pdf_interface.h"

namespace cpt_b0_analysis {

class RaisedCosinePDF : public PDFInterface {
public:
	RaisedCosinePDF() {};
	RaisedCosinePDF(const RaisedCosinePDF &other) = default;
	RaisedCosinePDF(RaisedCosinePDF &&) = default;
	~RaisedCosinePDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return IntCos; };
	int getComponentCount() { return 0; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	double IntCos = 1.0;
};

class GaussPDF : public PDFInterface {
public:
	GaussPDF() {};
	GaussPDF(const GaussPDF &other) = default;
	GaussPDF(GaussPDF &&) = default;
	~GaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return IntGaus; };
	int getComponentCount() { return 0; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	double IntGaus = 1.0;
};

class DoubleGaussPDF : public PDFInterface {
public:
	DoubleGaussPDF() {};
	DoubleGaussPDF(const DoubleGaussPDF &other) = default;
	DoubleGaussPDF(DoubleGaussPDF &&) = default;
	~DoubleGaussPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return 0.0; };
	int getComponentCount() { return 2; };	// 2 x Gauss
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	GaussPDF Gauss1;
	GaussPDF Gauss2;
};

class SkewNormalPDF : public PDFInterface {
public:
	SkewNormalPDF() {};
	SkewNormalPDF(const SkewNormalPDF &other) = default;
	SkewNormalPDF(SkewNormalPDF &&) = default;
	~SkewNormalPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return IntSkewNorm; };
	int getComponentCount() { return 0; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	double IntSkewNorm = 1.0;
};

class CrystalBallPDF : public PDFInterface {
public:
	CrystalBallPDF() {};
	CrystalBallPDF(const CrystalBallPDF &other) = default;
	CrystalBallPDF(CrystalBallPDF &&) = default;
	~CrystalBallPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return IntCB; };
	int getComponentCount() { return 0; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	double IntCB = 1.0;
};

class JohnsonPDF : public PDFInterface {
public:
	JohnsonPDF() {};
	JohnsonPDF(const JohnsonPDF &other) = default;
	JohnsonPDF(JohnsonPDF &&) = default;
	~JohnsonPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return IntJSU; };
	int getComponentCount() { return 0; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	double IntJSU = 1.0;
};

class DoubleSidedCrystalBallPDF : public PDFInterface {
public:
	DoubleSidedCrystalBallPDF() {};
	DoubleSidedCrystalBallPDF(const DoubleSidedCrystalBallPDF &other) = default;
	DoubleSidedCrystalBallPDF(DoubleSidedCrystalBallPDF &&) = default;
	~DoubleSidedCrystalBallPDF() {};

	void CalcIntegral(const double *par, double min, double max);
	double getIntegral() { return IntDCB; };
	int getComponentCount() { return 0; };
	double EvalPDF(const double *xx, const double *par, const int component = -1);

private:
	double IntDCB = 1.0;
};
}  // namespace cpt_b0_analysis

#endif
// vim: tabstop=4 softtabstop=0 noexpandtab shiftwidth=4
