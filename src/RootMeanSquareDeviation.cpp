#include "RootMeanSquareDeviation.h"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "AuxMath.h"

using namespace std;

RootMeanSquareDeviation::RootMeanSquareDeviation(){}

RootMeanSquareDeviation::~RootMeanSquareDeviation(){}

/// LEMBRAR DE CORRIGIR OS ZEROS DO LUMPACVIEW
double RootMeanSquareDeviation::rmsd(string file1, string file2)
{
	vector<double> set1 = readPoint(file1);
	set1 = rotateToZ0(set1);
	set1 = rotateToPlane(set1);
	set1 = mirrorY(set1);
	vector<double> set2 = readPoint(file2);
//	set2 = rotateToZ0(set2);
//	set2 = rotateToPlane(set2);
	set2 = mirrorY(set2);
#ifdef _DEBUG
	printXyz(file1 + ".xyz", set1);
	printXyz(file2 + ".xyz", set2);
#endif

	double sum = 0.0e0;
	for (size_t i = 0; i < set1.size(); i++)
		sum += (set1[i] - set2[i])*(set1[i] - set2[i]);

	sum /= (double)set1.size();

	sum = sqrt(sum);

	return sum;
}

/*
operacoes no repulsion:
roda tudo, levando o 1o ponto pro eixo z.
coloca o vetor no eixo z e roda o 2o ponto para o plano xz com x>0.
checa se o 3o ponto tem y positivo.
se sim, termina.
se nao, reflete os pontos no plano xz.
*/


/*
formato dos arquivos:
file
[
n pontos
x1
x2
x3
...
y1
y2
y3
...
z1
z2
z3
...
]

os dois arquivos precisam estar alinhados. C1 = C1, Oxi1 = Oxi2 e etc.

*/


vector<double> RootMeanSquareDeviation::readPoint(string fName)
{
	ifstream fPoint_(fName.c_str());
	int nPoints;
	fPoint_ >> nPoints;
	double aux;
	vector<double> points;
	for (int i = 0; i < nPoints; i++)
	{
		fPoint_ >> aux;
		points.push_back(aux);
	}
	fPoint_.close();
	return points;
}

vector<double> RootMeanSquareDeviation::rotateToZ0(const vector<double> &point)
{
	int nPoints = point.size() / 3;
	double x = point[0];
	double y = point[0 + nPoints];
	double z = point[0 + 2 * nPoints];

	AuxMath auxMath_;

	double angle = auxMath_.angleFrom3Points(
		x, y, z,
		0.0e0, 0.0e0, 0.0e0,
		0.0e0, 0.0e0, 1.0e0);

	vector<double> normal = auxMath_.normalVectorFrom3Points(
		x, y, z,
		0.0e0, 0.0e0, 0.0e0,
		0.0e0, 0.0e0, 1.0e0);

	vector< vector<double> > rot = auxMath_.rotationMatrix(normal[0], normal[1], normal[2], -angle);

	vector<double> newPoints(3 * nPoints);
	vector<double> aux;
	for (int i = 0; i < nPoints; i++)
	{
		aux = auxMath_.matrixXVector(rot, point[i], point[i + nPoints], point[i + 2 * nPoints]);
		newPoints[i] = aux[0];
		newPoints[i + nPoints] = aux[1];
		newPoints[i + 2 * nPoints] = aux[2];
	}
	return newPoints;
}

vector<double> RootMeanSquareDeviation::rotateToPlane(const vector<double> &point)
{
	AuxMath auxMath_;
	int nPoints = point.size() / 3;
	int point2Chosen;
	//ele precisa estar mais proximo do zero.
	double rMin = 1.0e99;
	double r;
	for (size_t i = 1; i < nPoints; i++)
	{
		r = auxMath_.norm(
			point[i] - point[0],
			point[i + nPoints] - point[0 + nPoints],
			point[i + 2 * nPoints] - point[0 + 2 * nPoints]
			);
		cout << setprecision(16) << r << endl;
		if (r < rMin)
		{
			point2Chosen = i;
			rMin = r;
		}
	}

	double x = point[point2Chosen];
	double y = point[point2Chosen + nPoints];
	double z = point[point2Chosen + 2 * nPoints];

	double angle = auxMath_.angleFrom3Points(
		x, y, 0.0e0,
		0.0e0, 0.0e0, 0.0e0,
		1.0e0, 0.0e0, 0.0e0);

	vector<double> normal(3);
	normal[0] = 0.0e0;
	normal[1] = 0.0e0;
	normal[2] = 1.0e0;

	if ((x < 0.0e0) && (y > 0))
		angle *= -1;

	vector< vector<double> > rot = auxMath_.rotationMatrix(normal[0], normal[1], normal[2], angle);

	vector<double> newPoints(3 * nPoints);
	vector<double> aux;
	for (int i = 0; i < nPoints; i++)
	{
		aux = auxMath_.matrixXVector(rot, point[i], point[i + nPoints], point[i + 2 * nPoints]);
		newPoints[i] = aux[0];
		newPoints[i + nPoints] = aux[1];
		newPoints[i + 2 * nPoints] = aux[2];
	}
	return newPoints;
}



vector<double> RootMeanSquareDeviation::mirrorY(const vector<double> &point)
{
	int nPoints = point.size() / 3;
	vector<double> newPoints(3 * nPoints);

	if (point[2 + nPoints] < 0.0e0)
	{
		for (int i = 0; i < nPoints; i++)
		{
			newPoints[i] = point[i];
			newPoints[i + nPoints] = -point[i + nPoints];
			newPoints[i + 2 * nPoints] = point[i + 2 * nPoints];
		}
	}
	else
		newPoints = point;

	return newPoints;
}

void RootMeanSquareDeviation::printXyz(string fName, const vector<double> &points)
{
	int nPoints = points.size() / 3;
	ofstream xyz_(fName.c_str());
	xyz_ << nPoints << endl;
	xyz_ << "t" << endl;
	double x, y, z;
	for (int i = 0; i < nPoints; i++)
	{
		x = points[i];
		y = points[i + nPoints];
		z = points[i + 2 * nPoints];
		//if (abs(x) < 1.0e-4)
			//x = 0.0e0;
		//if (abs(y) < 1.0e-4)
			//y = 0.0e0;
		//if (abs(z) < 1.0e-4)
			//z = 0.0e0;

		xyz_ << "H  " 
			<< setfill(' ')  << setw(20)	<< fixed << setprecision(16)  
			<< x << "  "
			<< setfill(' ') << setw(20) << fixed << setprecision(16)
			<< y << "  "
			<< setfill(' ') << setw(20) << fixed << setprecision(16)
			<< z << "  " << endl;
	}
}
