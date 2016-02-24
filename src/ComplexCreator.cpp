#include "ComplexCreator.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Ligand.h"
#include "Points.h"
#include "AuxMath.h"

using namespace std;

ComplexCreator::ComplexCreator(
	vector<Ligand> &allLigands_in,
	string metalName_in,
	string metalParams_in,
	string projectName_in
	)
	:allLigands(allLigands_in)
{
	metalName = metalName_in;
	metalParams = metalName_in;
	projectName = projectName_in;
}

ComplexCreator::~ComplexCreator(){}

bool ComplexCreator::start()
{
	int sumChelation = orderAllLigands();

	if (sumChelation > maxChelation)
		return false;

	vector<double> points = getPoints(sumChelation);

	stretchPoints(points);

	setInitialPosition(points);

#ifdef _DEBUG
	ofstream points_("points.xyz");
	int nPoints = (int)points.size() / 3;
	points_ << nPoints << endl << "t" << endl;
	for (int i = 0; i < nPoints; i++)
		points_ << "H  " << points[i] << "   "
			<< points[i + nPoints] << "   "
			<< points[i + 2 * nPoints] << endl;
	points_.close();

	printAllAtoms();
#endif 

	return true;
}


vector<double> ComplexCreator::getPoints(int totalChelation)
{
	Points allPoints;
	// 2 points - 0 - 5; 3 points - 6 - 15; 4 points - 16 - 28 ...
	int k = 0;
	for (int i = 3; i <= totalChelation; i++)
		k += (3 * (i - 1));

	size_t size = totalChelation;
	vector<double> out(3 * size);
	for (size_t i = 0; i < size; i++)
	{
		out[i] = allPoints.p[k + i];
		out[i + size] = allPoints.p[k + i + size];
		out[i + 2 * size] = allPoints.p[k + i + 2 * size];
	}
	return out;
}

vector<double> ComplexCreator::arrayToVector(
	const double * array_in, 
	size_t size)
{
	vector<double> out(3 * size);
	for (size_t i = 0; i < size; i++)
	{
		out[i] = array_in[i];
		out[i] = array_in[i + size];
		out[i] = array_in[i + 2 * size];
	}

	return out;
}


int ComplexCreator::orderAllLigands()
{
	size_t size = allLigands.size();
	vector<int> monoPositions;
	vector<int> biPositions;
	vector<int> triPositions;
	int chelation;
	int sumChelation = 0;
	for (size_t i = 0; i < size; i++)
	{
		chelation = allLigands[i].getChelation();
		sumChelation += chelation;
		switch (chelation)
		{
		case 1:
			monoPositions.push_back(i);
			break;
		case 2:
			biPositions.push_back(i);
			break;
		default:
			triPositions.push_back(i);
			break;
		}
	}

	vector<Ligand> tempLigand(size);
	int k = 0;
	for (size_t i = 0; i < triPositions.size(); i++)
	{
		tempLigand[k] = allLigands[triPositions[i]];
		k++;
	}
	for (size_t i = 0; i < biPositions.size(); i++)
	{
		tempLigand[k] = allLigands[biPositions[i]];
		k++;
	}
	for (size_t i = 0; i < monoPositions.size(); i++)
	{
		tempLigand[k] = allLigands[monoPositions[i]];
		k++;
	}

	allLigands = tempLigand;

	return sumChelation;
}


void ComplexCreator::setInitialPosition(
	const vector<double> &points)
{
	size_t nPoints = points.size() / 3;
	vector<bool> pointsTaken(nPoints);
	for (size_t k = 0; k < nPoints; k++)
		pointsTaken[k] = false;

	for (size_t i = 0; i < allLigands.size(); i++)
	{
		vector<double> translate = findGoodPoint(
			allLigands[i].getChelation(),
			points,
			pointsTaken);

		allLigands[i].translateLigand(translate[0], translate[1], translate[2]);
		
		allLigands[i].rotateToCenter();
	}
}


vector<double> ComplexCreator::findGoodPoint(
	int chelation,
	const vector<double> &points,
	vector<bool>& pointsTaken)
{
	size_t nPoints = points.size() / 3;
	int p1;
	for (size_t k = 0; k < nPoints; k++)
	{
		if (!pointsTaken[k])
		{
			p1 = k;
			break;
		}
	}
	pointsTaken[p1] = true;
	vector<int> pLig(1);
	pLig[0] = p1;

	// p2: closest - p3: closest to mean point.
	if (chelation > 1)
	{
		pointsTaken[p1] = true;
		int p2 = closestPoint(
			points[p1],
			points[p1 + nPoints],
			points[p1 + 2 * nPoints],
			points,
			pointsTaken);
		pLig.push_back(p2);
		if (chelation > 2)
		{
			double xm = (points[p1] + points[p2]) / 2.0e0;
			double ym = (points[p1 + nPoints] + points[p2 + nPoints]) / 2.0e0;
			double zm = (points[p1 + 2 * nPoints] + points[p2 + 2 * nPoints]) / 2.0e0;
			int p3 = closestPoint(
				xm, ym, zm,
				points,
				pointsTaken);
			pLig.push_back(p3);
		}
	}

	vector<double> ligandPoint(3);
	ligandPoint[0] = points[pLig[0]];
	ligandPoint[1] = points[pLig[0] + nPoints];
	ligandPoint[2] = points[pLig[0] + 2 * nPoints];
	// chelation 2: mean - chelation 3: barycenter
	if (chelation == 2)
	{
		ligandPoint[0] = (ligandPoint[0] + points[pLig[1]]) / 2.0e0;
		ligandPoint[1] = (ligandPoint[1] + points[pLig[1] + nPoints]) / 2.0e0;
		ligandPoint[2] = (ligandPoint[2] + points[pLig[1] + 2 * nPoints]) / 2.0e0;
	}
	else if (chelation == 3)
	{
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

	return ligandPoint;
}


int ComplexCreator::closestPoint(
	double x, double y, double z,
	const vector<double> &points,
	std::vector<bool>& pointsTaken)
{
	AuxMath auxMath_;
	size_t nPoints = points.size() / 3;
	double closest = 100.0e0;
	int iClose = 0;
	double r;
	for (size_t i = 0; i < nPoints; i++)
	{
		r = auxMath_.norm(
			x - points[i],
			y - points[i + nPoints],
			z - points[i + 2 * nPoints]);

		if ((!pointsTaken[i]) &&
			(r < closest))
		{
			closest = r;
			iClose = i;
		}
	}
	pointsTaken[iClose] = true;
	return iClose;
}


void ComplexCreator::stretchPoints(vector<double> &points)
{
	size_t size = points.size();
	for (size_t i = 0; i < size; i++)
		points[i] *= stretchDistance;
}


void ComplexCreator::printAllAtoms()
{
	int allAtoms = 0;
	size_t size = allLigands.size();
	for (size_t k = 0; k < size; k++)
		allAtoms += allLigands[k].getNatoms();

	ofstream xyzAll("xyzAll.xyz");
	xyzAll << allAtoms + 1 << endl
		<< projectName << endl;
	xyzAll << metalName << "   0.00000    0.00000     0.00000" << endl;
	for (size_t i = 0; i < size; i++)
		allLigands[i].printLigand(xyzAll);

}



/*
int nPoints = size / 3;
vector< vector<double> > distanceMatrix(nPoints);
for (size_t l = 0; l < nPoints; l++)
	distanceMatrix[l].resize(nPoints);

for (size_t i = 0; i < (nPoints - 1); i++)
{
	for (size_t j = i + 1; j < nPoints; j++)
	{
		distanceMatrix[i][j] = auxMath_.norm(
			points[i] - points[j],
			points[i + nPoints] - points[j + nPoints],
			points[i + 2 * nPoints] - points[j + 2 * nPoints]
			);
		distanceMatrix[j][i] = distanceMatrix[j][i];
	}
}
*/
