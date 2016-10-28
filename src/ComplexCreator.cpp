#include "ComplexCreator.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Ligand.h"
#include "Points.h"
#include "AuxMath.h"
#include "Fitness.h"

using namespace std;

ComplexCreator::ComplexCreator(vector<Ligand> &allLigands_in)
	:allLigands(allLigands_in)
{
}

ComplexCreator::ComplexCreator(
	vector<Ligand> &allLigands_in,
	int maxChelation_in,
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
	maxAlfaAngle = maxAlfaAngle_in;
	maxBetaAngle = maxBetaAngle_in;
	saMaxIterations = saMaxIterations_in;
	saTemperatureUpdate = saTemperatureUpdate_in;
	saInitialTemperature = saInitialTemperature_in;
	saAcceptance = saAcceptance_in;
}

ComplexCreator::~ComplexCreator(){}

void ComplexCreator::calculateAllAngles(int nPoints)
{
	AuxMath auxMath_;
	vector<double> points = getPoints(nPoints);
	double xi, yi, zi, xj, yj, zj;
	ofstream of_("todosangulos.csv");
	for (int i = 0; i < nPoints - 1; i++)
	{
		for (int j = i + 1; j < nPoints; j++)
		{
			xi = points[i];
			yi = points[i + nPoints];
			zi = points[i + 2 * nPoints];

			xj = points[j];
			yj = points[j + nPoints];
			zj = points[j + 2 * nPoints];

			double angle = auxMath_.angleFrom3Points(xi, yi, zi, 0.0e0, 0.0e0, 0.0e0, xj, yj, zj);
			of_ << "i: " << i << " j: " << j << " angle:     ;  " << angle*(180/auxMath_._pi) << endl;
		}
	}
	of_.close();
}

bool ComplexCreator::start()
{
	vector<int> ligandsPermutation;
	return start(ligandsPermutation);
}

bool ComplexCreator::start(vector<int> & ligandsPermutation)
{
	int sumChelation = orderAllLigands();

	if (sumChelation > maxChelation)
		return false;

	vector<double> points = getPoints(sumChelation);

	setInitialPosition(points, ligandsPermutation); //view fredmudar

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
	const vector<double> &points,
	vector<int> & ligandsPermutation)
{
	size_t nPoints = points.size() / 3;
	vector<bool> pointsTaken(nPoints);
	for (size_t k = 0; k < nPoints; k++)
		pointsTaken[k] = false;

	if (ligandsPermutation.size() == 0)
	{
		for (size_t i = 0; i < allLigands.size(); i++)
		{
			vector<int> pointsOverLigand = findGoodPoint(
				allLigands[i].getChelation(),
				points,
				pointsTaken);

			allLigands[i].placeLigandOnPoins(
				pointsOverLigand,
				points);

			allLigands[i].rotateOverReferencePoints();
		}
	}
	else
	{
		int k = 0;
		for (size_t i = 0; i < allLigands.size(); i++)
		{
			vector<int> pointsOverLigand(allLigands[i].getChelation());

			for (size_t j = 0; j < allLigands[i].getChelation(); j++)
			{
				pointsOverLigand[j] = ligandsPermutation[k];
				k++;
			}

			allLigands[i].placeLigandOnPoins(
				pointsOverLigand,
				points);

			allLigands[i].rotateOverReferencePoints();
		}
	}
}


vector<int> ComplexCreator::findGoodPoint(
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

	return pLig;
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
	printAllAtoms("zinitialComplexBeforeAnnealing.xyz", x0);
	ofstream annealingInfo_("AnnealingInfo.log");
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
			prob = exp((f0 - f) / cTemp);
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
			finalI = i;
		}

#ifdef _DEBUG
		annealingInfo_ << "i:  " << i << "  fMin:  "
			<< setprecision(16) 
			<< fMin
			<< "  cTemp:  " << cTemp << endl;
#endif
	}

#ifdef _DEBUG
	printAllAtoms("zcomplexAfterAnnealing.xyz", xMin);
	annealingInfo_.close();
#endif

	vector<CoordXYZ> allAtoms;
	for (size_t i = 0; i < xMin.size(); i++)
	{
		vector<CoordXYZ> atomsLigand = xMin[i].getAllAtoms();
		allAtoms.insert(allAtoms.end(), atomsLigand.begin(), atomsLigand.end());
	}
	return allAtoms;
}

std::vector<Ligand> ComplexCreator::getLigandsCreated() const
{
	return allLigands;
}


double ComplexCreator::calculateAllFit(vector<Ligand> & ligands)
{
	double allFit = 0.0e0;
	Fitness fit_;
	for (size_t i = 0; i < (allLigands.size() - 1); i++)
	{
		for (size_t j = i + 1; j < allLigands.size(); j++	)
		{
			vector<CoordXYZ> atomsI = ligands[i].getAllAtoms();
			vector<CoordXYZ> atomsJ = ligands[j].getAllAtoms();
			allFit += fit_.calculateFit(atomsI,atomsJ);
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
		randomRot[0] = auxMath_.fRand(0.001, 1.0e0);
		randomRot[1] = auxMath_.fRand(0.001, 1.0e0);
		randomRot[2] = auxMath_.fRand(0.001, 1.0e0);
		auxMath_.normalize(randomRot);
		alfa = auxMath_.fRand(0, maxAlfaAngle);

		ligands[i].genericRotation(randomRot[0], randomRot[1], randomRot[2], alfa);

		beta = auxMath_.fRand(0, maxBetaAngle);
		ligands[i].rotateOnX2(beta);
	}
}



void ComplexCreator::printAllAtoms(std::string xyzName, vector<Ligand> & ligands)
{
	int allAtoms = 0;
	size_t size = ligands.size();
	for (size_t k = 0; k < size; k++)
		allAtoms += ligands[k].getNatoms();

	ofstream xyzAll(xyzName.c_str());
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

