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

	//SORTING LIST
	// calculating spherical teta and fi
	AuxMath auxMath_;
	int nAtoms = set1.size() / 3;
	vector<double> teta(nAtoms);
	vector<double> fi(nAtoms);
	float x, y, z;
	for (int i = 1; i < nAtoms; i++)
	{
		x = set1[i];
		y = set1[i + nAtoms];
		z = set1[i + 2 * nAtoms];
		if ((abs(x) < 1.0e-10) && (abs(y) < 1.0e-10) && (z < 0))
		{
			teta[i] = 2 * auxMath_._pi;
			fi[i] = 0.0e0;
		}
		else
		{
			fi[i] = atan2(y , x);
			if (fi[i] < 0)
				fi[i] += 2.0e0 * auxMath_._pi;

			teta[i] = (float)acos(z / sqrt(x*x + y*y + z*z));
		}
	}

	// sorting
	// fi crescente
	// teta crescente
	vector<double> auxPoints = set1;
	double auxX, auxY, auxZ;
	double aux;
	double tetaMin, fiMin;
	int jMin;
	for (int i = 1; i < (nAtoms - 1); i++)
	{
		tetaMin = teta[i];
		fiMin = fi[i];
		jMin = i;
		for (int j = i + 1; j < nAtoms; j++)
		{
			if (abs(teta[j] - tetaMin) < 1.0e-6)
			{
				if (fi[j] < fiMin)
				{
					tetaMin = teta[j];
					fiMin = fi[j];
					jMin = j;
				}
			}
			else if(teta[j] < tetaMin)
			{
				tetaMin = teta[j];
				fiMin = fi[j];
				jMin = j;
			}
		}
		// substitui i em j --- j em i.
		aux = teta[jMin];
		teta[jMin] = teta[i];
		teta[i] = aux;

		aux = fi[jMin];
		fi[jMin] = fi[i];
		fi[i] = aux;

		auxX = auxPoints[jMin];
		auxY = auxPoints[jMin + nAtoms];
		auxZ = auxPoints[jMin + 2 * nAtoms];
		auxPoints[jMin] = auxPoints[i];
		auxPoints[jMin + nAtoms] = auxPoints[i + nAtoms];
		auxPoints[jMin + 2 * nAtoms] = auxPoints[i + 2 * nAtoms];
		auxPoints[i] = auxX;
		auxPoints[i + nAtoms] = auxY;
		auxPoints[i + 2 * nAtoms] = auxZ;
	}

	double sum = 0.0e0;

/*
SEGUNDA MOLECULA
	vector<double> set2 = readPoint(file2);
	set2 = rotateToZ0(set2);
	set2 = rotateToPlane(set2);
	set2 = mirrorY(set2);
	for (size_t i = 0; i < set1.size(); i++)
		sum += (set1[i] - set2[i])*(set1[i] - set2[i]);
	sum /= (double)set1.size();
	sum = sqrt(sum);
*/
	  
#ifdef _DEBUG
	printXyz(file1 + ".xyz", auxPoints);
//	printXyz(file2 + ".xyz", set2);
#endif

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
	vector<double> points(nPoints);


	int format = 1;


	if (format == 0)
	{
		for (int i = 0; i < nPoints; i++)
		{
			fPoint_ >> aux;
			points[i] = aux;
		}
	}
	else
	{
		int natoms = nPoints / 3;
		double aux2, aux3;
		for (int i = 0; i < natoms; i++)
		{
			fPoint_ >> aux;
			fPoint_ >> aux2;
			fPoint_ >> aux3;
			points[i] = aux;
			points[i + natoms] = aux2;
			points[i + 2 * natoms] = aux3;
		}
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

	if ((abs(x) < 1.0e-12) && (abs(y) < 1.0e-12))
		return point;

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
	for (size_t i = 1; (int)i < nPoints; i++)
	{
		r = auxMath_.norm(
			point[i] - point[0],
			point[i + nPoints] - point[0 + nPoints],
			point[i + 2 * nPoints] - point[0 + 2 * nPoints]
			);
		cout << "r:  " << setprecision(16) << r << endl;
		if (r < rMin)
		{
			point2Chosen = i;
			rMin = r;
		}
	}
//	cin.get();

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
	vector<double> pointsReordering = newPoints;
	pointsReordering[1] = newPoints[point2Chosen];
	pointsReordering[1 + nPoints] = newPoints[point2Chosen + nPoints];
	pointsReordering[1 + 2 * nPoints] = newPoints[point2Chosen + 2 * nPoints];
	pointsReordering[point2Chosen] = newPoints[1];
	pointsReordering[point2Chosen + nPoints] = newPoints[1 + nPoints];
	pointsReordering[point2Chosen + 2 * nPoints] = newPoints[1 + 2 * nPoints];

	return pointsReordering;
}



vector<double> RootMeanSquareDeviation::mirrorY(const vector<double> &point)
{
	AuxMath auxMath_;
	int nPoints = point.size() / 3;
	int point3Chosen;
	//maior Y
	double yMax = -1.0e99;
	double y;
	for (size_t i = 2; (int)i < nPoints; i++)
	{
		y = point[i + nPoints];
		if (y > yMax)
		{
			point3Chosen = i;
			yMax = y;
		}
	}
	
	vector<double> newPoints = point;

	if (point[point3Chosen + nPoints] < 0)
	{		
		for (int i = 0; i < nPoints; i++)
		{
			newPoints[i] = point[i];
			newPoints[i + nPoints] = -point[i + nPoints];
			newPoints[i + 2 * nPoints] = point[i + 2 * nPoints];
		}
	}

	vector<double> pointsReordering = newPoints;
	pointsReordering[2] = newPoints[point3Chosen];
	pointsReordering[2 + nPoints] = newPoints[point3Chosen + nPoints];
	pointsReordering[2 + 2 * nPoints] = newPoints[point3Chosen + 2 * nPoints];
	pointsReordering[point3Chosen] = newPoints[2];
	pointsReordering[point3Chosen + nPoints] = newPoints[2 + nPoints];
	pointsReordering[point3Chosen + 2 * nPoints] = newPoints[2 + 2 * nPoints];

	return pointsReordering;
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
			<< setfill(' ')  << setw(15)	<< fixed << setprecision(8)  
			<< x << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< y << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< z << "  " << endl;
	}
}





/*			if (abs(x) < 1.0e-10)
{
if(y > 0)
fi[i] = auxMath_._pi / 2.0e0;
else
fi[i] = (3.0e0 * auxMath_._pi) / 2.0e0;
}
else
*/
