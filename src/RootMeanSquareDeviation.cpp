#include "RootMeanSquareDeviation.h"

#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "AuxMath.h"
#include "Coordstructs.h"

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
		x = (float)set1[i];
		y = (float)set1[i + nAtoms];
		z = (float)set1[i + 2 * nAtoms];
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
	  
#ifdef _DEBUG
	printXyz(file1 + ".xyz", auxPoints);
#endif

	return sum;
}


double RootMeanSquareDeviation::rmsOverlay(string molName1, string molName2)
{
	vector<CoordXYZ> mol1 = readCoord(molName1);

	vector<CoordXYZ> mol2 = readCoord(molName2);

	return rmsOverlay(mol1, mol2);
}


double RootMeanSquareDeviation::rmsOverlay(vector<CoordXYZ> & mol1, vector<CoordXYZ> & mol2)
{
	AuxMath auxMath_;
	int nAtoms = mol1.size();
	if (mol1.size() != mol2.size())
	{
		cout << "error on RootMeanSquareDeviation - different molecules" << endl;
		exit(1);
	}

	//translate to center of mass
	CoordXYZ cm1, cm2;
	cm1.x = 0.0e0;
	cm1.y = 0.0e0;
	cm1.z = 0.0e0;
	cm2.x = 0.0e0;
	cm2.y = 0.0e0;
	cm2.z = 0.0e0;
	for (int i = 0; i < nAtoms; i++)
	{
		cm1.x += mol1[i].x;
		cm1.y += mol1[i].y;
		cm1.z += mol1[i].z;
		cm2.x += mol2[i].x;
		cm2.y += mol2[i].y;
		cm2.z += mol2[i].z;
	}
	cm1.x /= (double)nAtoms;
	cm1.y /= (double)nAtoms;
	cm1.z /= (double)nAtoms;
	cm2.x /= (double)nAtoms;
	cm2.y /= (double)nAtoms;
	cm2.z /= (double)nAtoms;
	for (int i = 0; i < nAtoms; i++)
	{
		mol1[i].x -= cm1.x;
		mol1[i].y -= cm1.y;
		mol1[i].z -= cm1.z;
	  	mol2[i].x -= cm2.x;
		mol2[i].y -= cm2.y;
		mol2[i].z -= cm2.z;
	}

#ifdef _DEBUG
	printXyz("mol1-at-cm.xyz", mol1);
	printXyz("mol2-at-cm.xyz", mol2);
#endif

	// agora eu quero rodar essas estruturas até que elas se superponham.
	// posso varrer os angulos teta e fi. Isso me da vetores apontando
	// para todas as direcoes.
	// Pra cada um desses vetores eu posso girar a vontade um psi.
	// outra alternativa é ajustar os primeiros atomos e deus salve a america.

	// a ideia mais pratica e pegar alguns atomos, rodar ele pra la e fazer algumas
	// so rodar pros 3 pontos e pronto.
	
	vector<double> normalMol1, normalMol2;
	bool colinear = true;
	int atom2 = 1;
	while (colinear)
	{
		normalMol1 = auxMath_.normalVectorFrom3Points(
			mol1[0].x, mol1[0].y, mol1[0].z,
			0.0e0, 0.0e0, 0.0e0,
			mol1[atom2].x, mol1[atom2].y, mol1[atom2].z);

		normalMol2 = auxMath_.normalVectorFrom3Points(
			mol2[0].x, mol2[0].y, mol2[0].z,
			0.0e0, 0.0e0, 0.0e0,
			mol2[atom2].x, mol2[atom2].y, mol2[atom2].z);

		if (
			((abs(normalMol1[0]) < 1.0e-12) &&
				(abs(normalMol1[1]) < 1.0e-12) &&
				(abs(normalMol1[2]) < 1.0e-12)) ||
			((abs(normalMol2[0]) < 1.0e-12) &&
				(abs(normalMol2[1]) < 1.0e-12) &&
				(abs(normalMol2[2]) < 1.0e-12))
			)
		{
			atom2++;
			if (atom2 > (int)(mol1.size() - 1))
			{
				cout << "LINEAR MOLECULES ARE NOT IMPLEMENTED - CONTACT DEVELOPERS" << endl;
				exit(1);
			}
		}
		else
			break;
	}


	double angle = auxMath_.angleFrom3Points(
		normalMol1[0], normalMol1[1], normalMol1[2],
		0.0e0, 0.0e0, 0.0e0,
		normalMol2[0], normalMol2[1], normalMol2[2]);

	vector<double> referenceRotation = auxMath_.normalVectorFrom3Points(
		normalMol1[0], normalMol1[1], normalMol1[2],
		0.0e0, 0.0e0, 0.0e0,
		normalMol2[0], normalMol2[1], normalMol2[2]);

	if (!
		((abs(referenceRotation[0]) < 1.0e-12) &&
		(abs(referenceRotation[1]) < 1.0e-12) &&
		(abs(referenceRotation[2]) < 1.0e-12)))
	{
		rotateMol(mol2, referenceRotation[0], referenceRotation[1], referenceRotation[2], angle);

#ifdef _DEBUG
		printXyz("mol2-afterrotation-cm.xyz", mol2);
#endif
	}

	double angleFi = auxMath_.angleFrom3Points(
		mol1[0].x, mol1[0].y, mol1[0].z,
		0.0e0, 0.0e0, 0.0e0,
		mol2[0].x, mol2[0].y, mol2[0].z
		);

	rotateMol(mol2, normalMol1[0], normalMol1[1], normalMol1[2], angleFi);
	double auxRms = rms(mol1, mol2);

#ifdef _DEBUG
	printXyz("mol2-ultimate-rotation-cm.xyz", mol2);
#endif

	return auxRms;
}


vector<CoordXYZ> RootMeanSquareDeviation::readCoord(string fName)
{
	vector<CoordXYZ> mol;
	ifstream fCoord_(fName.c_str());
	int nAtoms;
	string dummy;
	fCoord_ >> nAtoms;
	fCoord_ >> dummy;
	mol.resize(nAtoms);
	for (int i = 0; i < nAtoms; i++)
	{
		fCoord_ >> mol[i].atomlabel
			>> mol[i].x
			>> mol[i].y
			>> mol[i].z;
	}
	fCoord_.close();
	return mol;
}

vector<double> RootMeanSquareDeviation::readPoint(string fName, int format)
{
	ifstream fPoint_(fName.c_str());
	int nPoints;
	fPoint_ >> nPoints;
	double aux;
	vector<double> points(nPoints);

	if (format == 0)
	{
		for (int i = 0; i < nPoints; i++)
		{
			fPoint_ >> aux;
			points[i] = aux;
		}
	}
	else if(format == 1)
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

		xyz_ << "H  " 
			<< setfill(' ')  << setw(15)	<< fixed << setprecision(8)  
			<< x << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< y << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< z << "  " << endl;
	}
}

void RootMeanSquareDeviation::printXyz(string fName, vector<CoordXYZ> & mol)
{
	int nAtoms = mol.size();
	ofstream xyz_(fName.c_str());
	xyz_ << nAtoms << endl;
	xyz_ << "t" << endl;
	for (int i = 0; i < nAtoms; i++)
	{
		xyz_ << mol[i].atomlabel << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< mol[i].x << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< mol[i].y << "  "
			<< setfill(' ') << setw(15) << fixed << setprecision(8)
			<< mol[i].z << "  " << endl;
	}
}

double RootMeanSquareDeviation::rms(vector<CoordXYZ>& mol1, vector<CoordXYZ>& mol2)
{
	double rms = 0.0e0;
	int nAtoms = mol1.size();
	for (int i = 0; i < nAtoms; i++)
	{
		rms += (mol1[i].x - mol2[i].x) * (mol1[i].x - mol2[i].x);
		rms += (mol1[i].y - mol2[i].y) * (mol1[i].y - mol2[i].y);
		rms += (mol1[i].z - mol2[i].z) * (mol1[i].z - mol2[i].z);
	}
	rms /= nAtoms;
	rms = sqrt(rms);
	return rms;
}

void RootMeanSquareDeviation::rotateMol(
	vector<CoordXYZ> &mol,
	double x,
	double y,
	double z,
	double angle)
{
	AuxMath auxMath_;
	vector< vector<double> > rotationMatrix = auxMath_.rotationMatrix(
		x,
		y,
		z,
		angle);
	for (size_t i = 0; i < mol.size(); i++)
	{
		vector<double> atomRotated = auxMath_.matrixXVector(
			rotationMatrix,
			mol[i].x,
			mol[i].y,
			mol[i].z);
		mol[i].x = atomRotated[0];
		mol[i].y = atomRotated[1];
		mol[i].z = atomRotated[2];
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


