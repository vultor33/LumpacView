#include "Ligand.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>

#include "Coordstructs.h"
#include "AuxMath.h"

using namespace std;

Ligand::Ligand() {}

Ligand::~Ligand() {}

void Ligand::setLigandCoordinates(
	vector<CoordXYZ>& coord_in,
	string titleInfo_in)
{
	coord = coord_in;
	titleInfo = titleInfo_in;
}


void Ligand::setLigandCoordinates(string ligandFileName)
{
	ifstream mol_(ligandFileName.c_str());
	if (!mol_.is_open())
	{
		cout << "ERRO ON Ligand::setLigandCoordinates - ligand file not found" << endl;
		exit(1);
	}
	string auxline;
	int nAtoms;
	getline(mol_, auxline);
	stringstream line0;
	line0 << auxline;
	line0 >> nAtoms;

	vector<CoordXYZ> coord;
	string titleInfo;
	coord.resize(nAtoms);

	int i = 0;
	getline(mol_, auxline);
	titleInfo = auxline;
	while (getline(mol_, auxline))
	{
		if (i == nAtoms)
			break;

		stringstream line;
		line << auxline;
		line >> coord[i].atomlabel
			>> coord[i].x
			>> coord[i].y
			>> coord[i].z;
		i++;
	}
	mol_.close();

	this->setLigandCoordinates(coord, titleInfo);
}



bool Ligand::initializeLigand()
{
	bool success1 = getInfoFromTitle();
	if (!success1)
		return false;

	bool sucess2;
	switch (chelation)
	{
	case 1:
		sucess2 = calculateMonodentate();
		break;

	case 2:
		sucess2 = calculateBidentate();
		break;

	case 3:
		sucess2 = calculateTridentate();
		break;

	default:
		sucess2 = false;
	}
	if (!sucess2)
		return false;

	// center on x1
	translateLigand(-X1.x, -X1.y, -X1.z);
	X1.x = 0.0e0;
	X1.y = 0.0e0;
	X1.z = 0.0e0;

	return sucess2;
}


bool Ligand::getInfoFromTitle()
{
	// can get other info like mass
	stringstream line;
	line << titleInfo;
	string chelationLetter;
	line >> formalName >> chelationLetter >> metalDistance;
	if (chelationLetter == "monodentate")
		chelation = 1;
	else if (chelationLetter == "bidentate")
		chelation = 2;
	else if (chelationLetter == "tridentate")
		chelation = 3;
	else
	{
		cout << "chelation not found - CHECK LIGAND FORMAT" << endl;
		cin.get();
		return false;
	}	
	return true;
}

bool Ligand::calculateMonodentate()
{
	if ((int)coord.size() < chelation)
		return false;

	X1 = coord[0];

	if (coord.size() == 1)
		return true;

	vector<double> centroid(3);
	centroid[0] = 0.0e0;
	centroid[1] = 0.0e0;
	centroid[2] = 0.0e0;
	for (size_t i = 0; i < coord.size(); i++)
	{
		centroid[0] += coord[i].x;
		centroid[1] += coord[i].y;
		centroid[2] += coord[i].z;
	}
	centroid[0] /= coord.size();
	centroid[1] /= coord.size();
	centroid[2] /= coord.size();

	//vector x1 positivo pra ele q eu vou
	vector<double> direction(3);
	direction[0] = X1.x - centroid[0];
	direction[1] = X1.y - centroid[1];
	direction[2] = X1.z - centroid[2];

	AuxMath auxMath_;
	auxMath_.normalize(direction);
	X2.x = direction[0];
	X2.y = direction[1];
	X2.z = direction[2];

#ifdef _DEBUG
	printXyzLigandDirection();
#endif
	return true;
}

bool Ligand::calculateBidentate()
{
	srand(3);
	if (coord.size() < 2)
	{
		cout << "Need at least two atoms" << endl;
		return false;
	}

	AuxMath auxMath_;
	CoordXYZ thirdReferenceAtom;
	X1.x = 0.5e0 * (coord[0].x + coord[1].x);
	X1.y = 0.5e0 * (coord[0].y + coord[1].y);
	X1.z = 0.5e0 * (coord[0].z + coord[1].z);
	if (coord.size() == 2)
	{
		thirdReferenceAtom.x = 0.0e0;
		thirdReferenceAtom.y = 0.0e0;
		thirdReferenceAtom.z = 1.0e0;
		bool notLi;
		int kLi = 0;
		do
		{
			kLi++;
			if (kLi > 1000000)
			{
				cout << "Infinity loop at Ligand::calculateBidentate - contact developers" << endl;
				return false;
			}
			thirdReferenceAtom.x += auxMath_.fRand(0.0, 0.5e0);
			thirdReferenceAtom.y += auxMath_.fRand(0.0, 0.5e0);
			thirdReferenceAtom.z += auxMath_.fRand(0.0, 0.5e0);
			double linearDependence1 = auxMath_.escalarProduct(
				-X1.x + coord[0].x,
				-X1.y + coord[0].y,
				-X1.z + coord[0].z,
				-X1.x + thirdReferenceAtom.x,
				-X1.y + thirdReferenceAtom.y,
				-X1.z + thirdReferenceAtom.z
				);
			notLi = (abs(linearDependence1) < 1.0e-6);
		} while (notLi);
	}
	else
	{
		thirdReferenceAtom = coord[2];
	}

	// rotate third atom point to 90 degrees
	vector<double> normal = auxMath_.normalVectorFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		X1.x, X1.y, X1.z,
		thirdReferenceAtom.x, thirdReferenceAtom.y, thirdReferenceAtom.z
		);
	double angle = auxMath_.angleFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		X1.x, X1.y, X1.z,
		thirdReferenceAtom.x, thirdReferenceAtom.y, thirdReferenceAtom.z
		);

	// quando e menos e quando e mais.
	double rotAngle = (-angle + auxMath_._pi / 2.0e0);

	vector< vector<double> > rot = auxMath_.rotationMatrix(normal[0], normal[1], normal[2], rotAngle);

	vector<double> direction = auxMath_.matrixXVector(
		rot, 
		-X1.x + thirdReferenceAtom.x,
		-X1.y + thirdReferenceAtom.y,
		-X1.z + thirdReferenceAtom.z);

	double angleEnd = auxMath_.angleFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		X1.x, X1.y, X1.z,
		X1.x + direction[0], X1.y + direction[1], X1.z + direction[2]
		);
	
	//angleEnd = auxMath_._pi / 2.0e0;
	if (abs(angleEnd - (auxMath_._pi / 2.0e0)) > 1.0e-4)
	{
		rotAngle *= -1.0e0; 
		rot = auxMath_.rotationMatrix(normal[0], normal[1], normal[2], rotAngle);
		direction = auxMath_.matrixXVector(
			rot,
			-X1.x + thirdReferenceAtom.x,
			-X1.y + thirdReferenceAtom.y,
			-X1.z + thirdReferenceAtom.z);
	}
	
	auxMath_.normalize(direction);
	X2.x = -direction[0];
	X2.y = -direction[1];
	X2.z = -direction[2];

#ifdef _DEBUG
printXyzLigandDirection();
#endif
return true;
}

bool Ligand::calculateTridentate()
{
	if (coord.size() < 3)
	{
		cout << "Need at least three atoms" << endl;
		return false;
	}

	AuxMath auxMath_;
	vector<double> geometricCenter = auxMath_.triangleCentroid(
		coord[0].x, coord[0].y, coord[0].z,
		coord[1].x, coord[1].y, coord[1].z,
		coord[2].x, coord[2].y, coord[2].z
		);
	X1.x = geometricCenter[0];
	X1.y = geometricCenter[1];
	X1.z = geometricCenter[2];
	vector<double> normal = auxMath_.normalVectorFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		coord[1].x, coord[1].y, coord[1].z,
		coord[2].x, coord[2].y, coord[2].z
		);
	X2.x = normal[0];
	X2.y = normal[1];
	X2.z = normal[2];

	if (coord.size() > 3)
	{
		vector<double> centroid(3);
		centroid[0] = 0.0e0;
		centroid[1] = 0.0e0;
		centroid[2] = 0.0e0;
		for (size_t i = 0; i < coord.size(); i++)
		{
			centroid[0] += coord[i].x;
			centroid[1] += coord[i].y;
			centroid[2] += coord[i].z;
		}
		centroid[0] /= coord.size();
		centroid[1] /= coord.size();
		centroid[2] /= coord.size();
		// if it points to centroid, change direction
		double angle = auxMath_.angleFrom3Points(
			X2.x + X1.x, X2.y + X1.y, X2.z + X1.z,
			X1.x, X1.y, X1.z,
			centroid[0], centroid[1], centroid[2]
			);
		if (angle < auxMath_._pi / 2.0e0)
		{
			X2.x *= -1.0e0;
			X2.y *= -1.0e0;
			X2.z *= -1.0e0;
		}
	}

#ifdef _DEBUG
	printXyzLigandDirection();
#endif

	return true;
}


void Ligand::translateLigand(double x, double y, double z)
{
	for (size_t i = 0; i < coord.size(); i++)
	{
		coord[i].x += x;
		coord[i].y += y;
		coord[i].z += z;
	}
	X1.x += x;
	X1.y += y;
	X1.z += z;
}

void Ligand::rotateToCenter()
{
	if (coord.size() == 1)
		return;

	// rotation about X1.
	AuxMath auxMath_;
	vector<double> normal = auxMath_.normalVectorFrom3Points(
		-X1.x, -X1.y, -X1.z,
		0.0e0, 0.0e0, 0.0e0,
		X2.x, X2.y, X2.z);

	if ((abs(normal[0]) < 1.0e-6) &&
		(abs(normal[1]) < 1.0e-6) &&
		(abs(normal[2]) < 1.0e-6))
		return;

	double angle = auxMath_.angleFrom3Points(
		-X1.x, -X1.y, -X1.z,
		0.0e0, 0.0e0, 0.0e0,
		X2.x, X2.y, X2.z);

	rotateOnX1(normal[0], normal[1], normal[2], -angle);
}

void Ligand::rotateOnX1(double vx, double vy, double vz, double ang)
{
	AuxMath auxMath_;
	vector< vector<double> > mrot = auxMath_.rotationMatrix(vx, vy, vz, ang);

	CoordXYZ oldX1 = X1;
	translateLigand(-X1.x, -X1.y, -X1.z);
	
	vector<double> auxRot;
	for (size_t i = 0; i < coord.size(); i++)
	{
		auxRot = auxMath_.matrixXVector(
			mrot, 
			coord[i].x, 
			coord[i].y, 
			coord[i].z);
		coord[i].x = auxRot[0];
		coord[i].y = auxRot[1];
		coord[i].z = auxRot[2];
	}
	auxRot = auxMath_.matrixXVector(
		mrot,
		X2.x,
		X2.y,
		X2.z);

	X2.x = auxRot[0];
	X2.y = auxRot[1];
	X2.z = auxRot[2];

	translateLigand(oldX1.x, oldX1.y, oldX1.z);
}

void Ligand::rotateOnX2(double beta)
{
	rotateOnX1(X2.x, X2.y, X2.z, beta);
}

void Ligand::genericRotation(double vx, double vy, double vz, double ang)
{
	AuxMath auxMath_;
	vector< vector<double> > mrot = auxMath_.rotationMatrix(vx, vy, vz, ang);
	vector<double> auxRot;
	for (size_t i = 0; i < coord.size(); i++)
	{
		auxRot = auxMath_.matrixXVector(
			mrot,
			coord[i].x,
			coord[i].y,
			coord[i].z);
		coord[i].x = auxRot[0];
		coord[i].y = auxRot[1];
		coord[i].z = auxRot[2];
	}
	auxRot = auxMath_.matrixXVector(
		mrot,
		X1.x,
		X1.y,
		X1.z);
	X1.x = auxRot[0];
	X1.y = auxRot[1];
	X1.z = auxRot[2];

	auxRot = auxMath_.matrixXVector(
		mrot,
		X2.x,
		X2.y,
		X2.z);
	X2.x = auxRot[0];
	X2.y = auxRot[1];
	X2.z = auxRot[2];
}

void Ligand::placeLigandOnPoins(vector<int>& pLig, const vector<double>& points_in)
{
	vector<double> points = points_in;
	vector<double> ligandPoint(3);
	int nPoints = points.size() / 3;	

	for (size_t i = 0; i < points.size(); i++)
		points[i] *= metalDistance;

	double x1, y1, z1;
	double x2, y2, z2;
	double x3, y3, z3;
	x1 = points[pLig[0]];
	y1 = points[pLig[0] + nPoints];
	z1 = points[pLig[0] + 2 * nPoints];
	sphereReferencePoints.push_back(x1);
	sphereReferencePoints.push_back(y1);
	sphereReferencePoints.push_back(z1);
	ligandPoint[0] = x1;
	ligandPoint[1] = y1;
	ligandPoint[2] = z1;
	// chelation 2: mean - chelation 3: barycenter
	if (chelation == 2)
	{
		x2 = points[pLig[1]];
		y2 = points[pLig[1] + nPoints];
		z2 = points[pLig[1] + 2 * nPoints];
		sphereReferencePoints.push_back(x2);
		sphereReferencePoints.push_back(y2);
		sphereReferencePoints.push_back(z2);
		ligandPoint[0] = (ligandPoint[0] + x2) / 2.0e0;
		ligandPoint[1] = (ligandPoint[1] + y2) / 2.0e0;
		ligandPoint[2] = (ligandPoint[2] + z2) / 2.0e0;
	}
	else if (chelation == 3)
	{
		x2 = points[pLig[1]];
		y2 = points[pLig[1] + nPoints];
		z2 = points[pLig[1] + 2 * nPoints];
		x3 = points[pLig[2]];
		y3 = points[pLig[2] + nPoints];
		z3 = points[pLig[2] + 2 * nPoints];
		sphereReferencePoints.push_back(x2);
		sphereReferencePoints.push_back(y2);
		sphereReferencePoints.push_back(z2);
		sphereReferencePoints.push_back(x3);
		sphereReferencePoints.push_back(y3);
		sphereReferencePoints.push_back(z3);
		ligandPoint[0] = (
			ligandPoint[0] +
			points[pLig[1]] +
			points[pLig[2]]) / 3.0e0;
		ligandPoint[1] = (
			ligandPoint[1] +
			points[pLig[1] + nPoints] +
			points[pLig[2] + nPoints]) / 3.0e0;
		ligandPoint[2] = (
			ligandPoint[2] +
			points[pLig[1] + 2 * nPoints] +
			points[pLig[2] + 2 * nPoints]) / 3.0e0;
	}

	translateLigand(ligandPoint[0], ligandPoint[1], ligandPoint[2]);

	rotateToCenter();

	if (chelation == 1)
		return;

	// rotating ligand to overlap points.
#ifdef _DEBUG
	ofstream placeLigand_("placeLigand0.xyz");
	placeLigand_ << coord.size() + chelation << endl << "t " << endl;
	printLigand(placeLigand_);
	placeLigand_ << "Au  " << x1 << "  " << y1 << "  " << z1 << endl;
	placeLigand_ << "Au  " << x2 << "  " << y2 << "  " << z2 << endl;
	if(chelation == 3)
		placeLigand_ << "Au  " << x3 << "  " << y3 << "  " << z3 << endl;

	placeLigand_.close();
#endif

	AuxMath auxMath_;
	int iDistanceMin = 0;
	double distanceMin, distance;

	distanceMin = auxMath_.norm(
		coord[0].x - x1,
		coord[0].y - y1,
		coord[0].z - z1);
	distanceMin += auxMath_.norm(
		coord[1].x - x2,
		coord[1].y - y2,
		coord[1].z - z2);
	if (chelation == 3)
	{
		distanceMin += auxMath_.norm(
			coord[2].x - x3,
			coord[2].y - y3,
			coord[2].z - z3);
	}
	for (int i = 1; i < 63; i++)
	{
		rotateOnX2(0.1e0);

		distance = auxMath_.norm(
			coord[0].x - x1,
			coord[0].y - y1,
			coord[0].z - z1);
		distance += auxMath_.norm(
			coord[1].x - x2,
			coord[1].y - y2,
			coord[1].z - z2);
		if (chelation == 3)
		{
			distance += auxMath_.norm(
				coord[2].x - x3,
				coord[2].y - y3,
				coord[2].z - z3);
		}
		if (distance < distanceMin)
		{
			distanceMin = distance;
			iDistanceMin = i;
		}
	}
	
	rotateOnX2(-(62-iDistanceMin) * 0.1);

#ifdef _DEBUG
	ofstream placeLigand2_;
	placeLigand2_.open("placeLigandEND.xyz");
	placeLigand2_ << coord.size() + chelation << endl << "t " << endl;
	printLigand(placeLigand2_);
	placeLigand2_ << "Au  " << x1 << "  " << y1 << "  " << z1 << endl;
	placeLigand2_ << "Au  " << x2 << "  " << y2 << "  " << z2 << endl;
	if (chelation == 3)
		placeLigand2_ << "Au  " << x3 << "  " << y3 << "  " << z3 << endl;
	placeLigand2_.close();
#endif

}

void Ligand::rotateOverReferencePoints(double angle)
{
	if (coord.size() == 1)
		return;

	AuxMath auxMath_;
	if (chelation == 1)
	{
		rotateOnX2(angle);
	}
	else if (chelation == 2)
	{
		rotateOnX2(auxMath_._pi);
	}
	else
	{
		rotateOnX2(auxMath_._pi / 2);
		double distance, distanceMin;
		int iDistanceMin = 0;
		distanceMin = auxMath_.norm(
			coord[0].x - sphereReferencePoints[0],
			coord[0].y - sphereReferencePoints[1],
			coord[0].z - sphereReferencePoints[2]);
		distanceMin += auxMath_.norm(
			coord[1].x - sphereReferencePoints[3],
			coord[1].y - sphereReferencePoints[4],
			coord[1].z - sphereReferencePoints[5]);
		distanceMin += auxMath_.norm(
			coord[2].x - sphereReferencePoints[6],
			coord[2].y - sphereReferencePoints[7],
			coord[2].z - sphereReferencePoints[8]);

		for (int i = 1; i < 10; i++)
		{
			rotateOnX2(0.1e0);
			distance = auxMath_.norm(
				coord[0].x - sphereReferencePoints[0],
				coord[0].y - sphereReferencePoints[1],
				coord[0].z - sphereReferencePoints[2]);
			distance += auxMath_.norm(
				coord[1].x - sphereReferencePoints[3],
				coord[1].y - sphereReferencePoints[4],
				coord[1].z - sphereReferencePoints[5]);
			distance += auxMath_.norm(
				coord[2].x - sphereReferencePoints[6],
				coord[2].y - sphereReferencePoints[7],
				coord[2].z - sphereReferencePoints[8]);

			if (distance < distanceMin)
			{
				distanceMin = distance;
				iDistanceMin = i;
			}
		}

		rotateOnX2(-(10 - iDistanceMin) * 0.1);
	}


#ifdef _DEBUG
	/*
	ofstream placeLigand_("placeLigand1.xyz");
	placeLigand_ << coord.size() + chelation << endl << "t " << endl;
	printLigand(placeLigand_);
	placeLigand_ << "Au  " << sphereReferencePoints[0] 
		<< "  " << sphereReferencePoints[1] 
		<< "  " << sphereReferencePoints[2] << endl;
	placeLigand_ << "Au  " << sphereReferencePoints[3] 
		<< "  " << sphereReferencePoints[4] 
		<< "  " << sphereReferencePoints[5] << endl;
	if (chelation == 3)
		placeLigand_ << "Au  " << sphereReferencePoints[6] 
		<< "  " << sphereReferencePoints[7] 
		<< "  " << sphereReferencePoints[8] << endl;

	placeLigand_.close();
	*/
#endif

}

double Ligand::distanceX1ToPoint(double x, double y, double z)
{
	return sqrt(
		(x - X1.x) * (x - X1.x)
		+ (y - X1.y) * (y - X1.y)
		+ (z - X1.z) * (z - X1.z));
}

vector<CoordXYZ> Ligand::getAllAtoms()
{
	return coord;
}

void Ligand::setNewCoordinates(vector<CoordXYZ>& newCoord)
{
	coord = newCoord;
	getInfoFromTitle();
	switch (chelation)
	{
	case 1:
		calculateMonodentate();
		break;

	case 2:
		calculateBidentate();
		break;

	case 3:
		calculateTridentate();
		break;
	}
}

void Ligand::printLigand(ofstream &out)
{
	for (size_t i = 0; i < coord.size(); i++)
	{
		out << coord[i].atomlabel << "   "
			<< setiosflags(ios::fixed) << setprecision(8) << coord[i].x << "   "
			<< setiosflags(ios::fixed) << setprecision(8) << coord[i].y << "   "
			<< setiosflags(ios::fixed) << setprecision(8) << coord[i].z << endl;
	}
}

void Ligand::printXyzLigandDirection(string inputName)
{
	vector<CoordXYZ> coordTemp = coord;
	CoordXYZ X1temp = X1;
	CoordXYZ X2temp = X2;
	X2temp.x += X1.x;
	X2temp.y += X1.y;
	X2temp.z += X1.z;
	X1temp.atomlabel = "Au";
	X2temp.atomlabel = "Cu";
	coordTemp.push_back(X1temp);
	coordTemp.push_back(X2temp);

	ofstream out_(inputName.c_str());
	size_t size = coordTemp.size();
	out_ << size << endl
		<< titleInfo << endl;
	for (size_t i = 0; i < size; i++)
	{
		out_ << coordTemp[i].atomlabel << "   "
			<< coordTemp[i].x << "   "
			<< coordTemp[i].y << "   "
			<< coordTemp[i].z << endl;
	}
	out_.close();
}


