#include "AuxMath.h"

#include <vector>

using namespace std;

vector<double> AuxMath::matrixXVector(vector< vector<double> > &M, double x, double y, double z)
{
	int size = M.size();
	vector<double> Mv(size);
	for (int i = 0; i < size; i++)
	{
		Mv[i] = M[i][0] * x + M[i][1] * y + M[i][2] * z;
	}
	return Mv;
}

double AuxMath::fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

double AuxMath::norm(const double x, const double y, const double z)
{
	return sqrt(x * x + y * y + z * z);
}

vector<double> AuxMath::vectorProduct(double x1, double y1, double z1, double x2, double y2, double z2)
{
	vector<double> prodvec(3);
	prodvec[0] = y1*z2 - y2*z1;
	prodvec[1] = z1*x2 - z2*x1;
	prodvec[2] = x1*y2 - x2*y1;

	return prodvec;
}

vector< vector<double> > AuxMath::rotationMatrix(double x, double y, double z, double teta)
{
	vector< vector<double> > rot(3);
	rot[0].resize(3);
	rot[1].resize(3);
	rot[2].resize(3);
	double ct = cos(teta);
	double st = sin(teta);
	rot[0][0] = x*x*(1.0e0 - ct) + ct;
	rot[0][1] = x*y*(1.0e0 - ct) - z*st;
	rot[0][2] = x*z*(1.0e0 - ct) + y*st;
	rot[1][0] = x*y*(1.0e0 - ct) + z*st;
	rot[1][1] = y*y*(1.0e0 - ct) + ct;
	rot[1][2] = y*z*(1.0e0 - ct) - x*st;
	rot[2][0] = x*z*(1.0e0 - ct) - y*st;
	rot[2][1] = y*z*(1.0e0 - ct) + x*st;
	rot[2][2] = z*z*(1.0e0 - ct) + ct;
	return rot;
}

std::vector<double> AuxMath::normalVectorFrom3Points(
	double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	double x3, double y3, double z3)
{
	double v1x, v1y, v1z, v2x, v2y, v2z;
	v1x = x1 - x2;
	v1y = y1 - y2;
	v1z = z1 - z2;
	v2x = x1 - x3;
	v2y = y1 - y3;
	v2z = z1 - z3;

	vector<double> n = vectorProduct(v1x, v1y, v1z, v2x, v2y, v2z);
	double module = norm(n[0], n[1], n[2]);
	n[0] /= module;
	n[1] /= module;
	n[2] /= module;

	return n;
}

std::vector<double> AuxMath::triangleCentroid(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
	vector<double> centroid(3);
	centroid[0] = 0.5e0 * (x1 + x2 + x3);
	centroid[1] = 0.5e0 * (y1 + y2 + y3);
	centroid[2] = 0.5e0 * (z1 + z2 + z3);

	return centroid;
}

double AuxMath::escalarProduct(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return x1*x2 + y1*y2 + z2*z2;
}


