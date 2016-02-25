#include "Ligand.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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

bool Ligand::initializeLigand()
{
	bool success1 = getInfoFromTitle();
	if (!success1)
		return false;

	// do math
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
	line >> chelationLetter;
	if (chelationLetter == "monodentate")
		chelation = 1;
	else if (chelationLetter == "bidentate")
		chelation = 2;
	else if (chelationLetter == "tridentate")
		chelation = 3;
	else
	{
		cout << "chelation not found" << endl;
		return false;
	}	
	return true;
}

bool Ligand::calculateMonodentate()
{
	if ((int)coord.size() < chelation)
		return false;

	X1 = coord[0];

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

	AuxMath auxCalc_;
	auxCalc_.normalize(direction);
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
	if (coord.size() < 3)
	{
		cout << "Need at least three atoms" << endl;
		return false;
	}

	X1.x = 0.5e0 * (coord[0].x + coord[1].x);
	X1.y = 0.5e0 * (coord[0].y + coord[1].y);
	X1.z = 0.5e0 * (coord[0].z + coord[1].z);

	AuxMath auxMath_;
	vector<double> normal = auxMath_.normalVectorFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		X1.x, X1.y, X1.z,
		coord[2].x, coord[2].y, coord[2].z
		);
	double angle = auxMath_.angleFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		X1.x, X1.y, X1.z,
		coord[2].x, coord[2].y, coord[2].z
		);

	vector< vector<double> > rot = auxMath_.rotationMatrix(normal[0], normal[1], normal[2], -angle + auxMath_._pi/2.0e0);

	vector<double> X2Points = auxMath_.matrixXVector(rot, coord[2].x, coord[2].y, coord[2].z);

	//vector x1 positivo pra ele q eu vou
	vector<double> direction(3);
	direction[0] = X1.x - X2Points[0];
	direction[1] = X1.y - X2Points[1];
	direction[2] = X1.z - X2Points[2];

	auxMath_.normalize(direction);
	X2.x = direction[0];
	X2.y = direction[1];
	X2.z = direction[2];

#ifdef _DEBUG
	printXyzLigandDirection();
#endif

	return true;
}

bool Ligand::calculateTridentate()
{
	if (coord.size() < 4)
	{
		cout << "Need at least four atoms" << endl;
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

	vector<double> normal = auxMath_.normalVectorFrom3Points(
		coord[0].x, coord[0].y, coord[0].z,
		coord[1].x, coord[1].y, coord[1].z,
		coord[2].x, coord[2].y, coord[2].z
		);
	X2.x = normal[0];
	X2.y = normal[1];
	X2.z = normal[2];

	// if it points to centroid, change direction
	double angle = auxMath_.angleFrom3Points(
		X2.x + X1.x, X2.y + X1.y, X2.z + X1.z,
		X1.x, X1.y, X1.z,
		centroid[0], centroid[1], centroid[2]
		);		

	if (angle < auxMath_._pi/2.0e0)
	{
		X2.x *= -1.0e0;
		X2.y *= -1.0e0;
		X2.z *= -1.0e0;
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
	// rotation about X1.
	AuxMath auxMath_;

	vector<double> normal = auxMath_.normalVectorFrom3Points(
		-X1.x, -X1.y, -X1.z,
		0.0e0, 0.0e0, 0.0e0,
		X2.x, X2.y, X2.z);

	double angle = auxMath_.angleFrom3Points(
		-X1.x, -X1.y, -X1.z,
		0.0e0, 0.0e0, 0.0e0,
		X2.x, X2.y, X2.z);

	rotateOnX1(normal[0], normal[1], normal[2], angle);
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

vector<CoordXYZ> Ligand::getAllAtoms()
{
	return coord;
}

void Ligand::printLigand(ofstream &out)
{
	for (size_t i = 0; i < coord.size(); i++)
	{
		out << coord[i].atomlabel << "   "
			<< coord[i].x << "   "
			<< coord[i].y << "   "
			<< coord[i].z << endl;
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


