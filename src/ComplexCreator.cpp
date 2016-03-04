#include "ComplexCreator.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Ligand.h"
#include "Points.h"
#include "AuxMath.h"
#include "Fitness.h"

using namespace std;

ComplexCreator::ComplexCreator(
	vector<Ligand> &allLigands_in,
	int maxChelation_in,
	double stretchDistance_in,
	int saMaxIterations_in,
	double maxAlfaAngle_in,
	double maxBetaAngle_in,
	double saTemperatureUpdate_in,
	double saInitialTemperature_in,
	double saAcceptance_in
	)
	:allLigands(allLigands_in)
{
	maxChelation = maxChelation_in;
	stretchDistance = stretchDistance_in;
	maxAlfaAngle = maxAlfaAngle_in;
	maxBetaAngle = maxBetaAngle_in;
	saMaxIterations = saMaxIterations_in;
	saTemperatureUpdate = saTemperatureUpdate_in;
	saInitialTemperature = saInitialTemperature_in;
	saAcceptance = saAcceptance_in;
}

ComplexCreator::~ComplexCreator(){}

bool ComplexCreator::start()
{
	int sumChelation = orderAllLigands();

	if (sumChelation > maxChelation)
		return false;

	vector<double> points = getPoints(sumChelation);

	stretchPoints(points);

	setInitialPosition(points); //view

#ifdef _DEBUG
	ofstream points_("points.xyz");
	int nPoints = (int)points.size() / 3;
	points_ << nPoints << endl << "t" << endl;
	for (int i = 0; i < nPoints; i++)
		points_ << "H  " << points[i] << "   "
			<< points[i + nPoints] << "   "
			<< points[i + 2 * nPoints] << endl;
	points_.close();
	printAllAtoms(allLigands);
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



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//////////////////  OPTIMIZING /////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

vector<CoordXYZ> ComplexCreator::simulatedAnnealing()
{
	AuxMath auxMath_;
	vector<Ligand> x0 = allLigands;
	double f0 = calculateAllFit(x0);
#ifdef _DEBUG
	printAllAtoms(x0);
	ofstream test_("teste.log");
#endif
	
	//Variable temperature - always 50/100
	double cTemp = f0 / saInitialTemperature; 
	int dices = 0;
	int accep = 0;

	vector<Ligand> xMin = x0;
	double fMin = f0;
	vector<Ligand> x;
	double f;
	double prob;
	double rand;
	for (int i = 0; i < saMaxIterations; i++)
	{
		x = x0;
		perturbOperations(x);
		f = calculateAllFit(x);
		if (f < f0)
		{
			x0 = x;
			f0 = f;
		}
		else
		{
			dices++;
			prob = exp((-f0 + f) / cTemp);
			rand = auxMath_.fRand(0, 1.0e0);
			if (prob > rand)
			{
				x0 = x;
				f0 = f;
				accep++;
			}
			if (((double)accep / (double)dices) > saAcceptance)
				cTemp -= saTemperatureUpdate * cTemp;
			else
				cTemp += saTemperatureUpdate * cTemp;
		}
		if (f < fMin)
		{
			xMin = x;
			fMin = f;
		}

#ifdef _DEBUG
		//printAllAtoms(xMin);
		test_ << "i:  " << i << "  fMin:  " << fMin
			<< "  cTemp:  " << cTemp << endl;
#endif
	}

#ifdef _DEBUG
	printAllAtoms(xMin);
	test_.close();
#endif

	vector<CoordXYZ> allAtoms;
	for (size_t i = 0; i < xMin.size(); i++)
	{
		vector<CoordXYZ> atomsLigand = xMin[i].getAllAtoms();
		allAtoms.insert(allAtoms.end(), atomsLigand.begin(), atomsLigand.end());
	}
	return allAtoms;
}


double ComplexCreator::calculateAllFit(vector<Ligand> & ligands)
{
	double allFit = 0.0e0;
	Fitness fit_;
	for (size_t i = 0; i < (allLigands.size() - 1); i++)
	{
		for (size_t j = i + 1; j < allLigands.size(); j++)
		{
			allFit += fit_.calculateFit(
				ligands[i].getAllAtoms(),
				ligands[j].getAllAtoms());
		}
	}
	return allFit;
}

void ComplexCreator::perturbOperations(vector<Ligand> & ligands)
{
	AuxMath auxMath_;
	vector<double> randomRot(3);
	double alfa, beta;
	for (size_t i = 0; i < ligands.size(); i++)
	{
		randomRot[0] = auxMath_.fRand(0, 1.0e0);
		randomRot[1] = auxMath_.fRand(0, 1.0e0);
		randomRot[2] = auxMath_.fRand(0, 1.0e0);
		auxMath_.normalize(randomRot);
		alfa = auxMath_.fRand(0, maxAlfaAngle);

		ligands[i].genericRotation(randomRot[0], randomRot[1], randomRot[2], alfa);

		beta = auxMath_.fRand(0, maxBetaAngle);
		ligands[i].rotateOnX2(beta);
	}
}



void ComplexCreator::printAllAtoms(vector<Ligand> & ligands)
{
	int allAtoms = 0;
	size_t size = ligands.size();
	for (size_t k = 0; k < size; k++)
		allAtoms += ligands[k].getNatoms();

	ofstream xyzAll("xyzAll.xyz");
	xyzAll << allAtoms + 1 << endl
		<< "useful title" << endl;
	xyzAll << "Eu  " << "   0.00000    0.00000     0.00000" << endl;
	for (size_t i = 0; i < size; i++)
		ligands[i].printLigand(xyzAll);

}






/*
TESTANDO SIMULATED ANNEALING
AuxMath auxMath_;

int maxIterations = 100;
double x0 = 3;
double f0 = x0 * x0;
vector<Ligand> ligand0 = allLigands;
double cTemp = f0; //Variable parameter is better

double xMin = x0;
double fMin = f0;
double x;
double f;
double prob;
for (int i = 0; i < maxIterations; i++)
{
x = x0 + 1.0e0 - auxMath_.fRand(0.0e0, 2.0e0);
f = x * x;

if (f < f0)
{
x0 = x;
f0 = f;
}
else
{
prob = exp((fMin - f) / cTemp);
if (prob > auxMath_.fRand(0, 1.0e0))
{
x0 = x;
f0 = f;
}
}
if (f < fMin)
{
xMin = x;
fMin = f;
}
cout << "x0=  " << x0 << endl;
}
cout << "MIN=  " << xMin << endl;
*/

