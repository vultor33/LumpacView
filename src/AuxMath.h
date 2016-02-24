#ifndef AUXMATH_H
#define AUXMATH_H

#include <vector>

class AuxMath
{
public:
	const double _pi = 3.14159265358979323846;

	std::vector<double> matrixXVector(std::vector< std::vector<double> > &M, double x, double y, double z);
	double fRand(double fMin, double fMax);
	double norm(const double x, const double y, const double z);
	void normalize(std::vector<double> &vector);
	std::vector<double> vectorProduct(double x1, double y1, double z1, double x2, double y2, double z2);
	std::vector< std::vector<double> > rotationMatrix(double x, double y, double z, double teta);
	std::vector<double> normalVectorFrom3Points(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
	std::vector<double> triangleCentroid(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
	double angleFrom3Points(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);

	double escalarProduct(double x1, double y1, double z1, double x2, double y2, double z2);

};

#endif

