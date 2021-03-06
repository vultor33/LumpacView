#include "ComplexCreator.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "Ligand.h"
#include "AuxMath.h"
#include "Fitness.h"
#include "CauchyIndex.h"

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

void ComplexCreator::calculateAllAngles(vector<CoordXYZ> &points)
{
	AuxMath auxMath_;
	double xi, yi, zi, xj, yj, zj;
	ofstream of_("todosangulos.csv");
	for (int i = 0; i < points.size() - 1; i++)
	{
		for (int j = i + 1; j < points.size(); j++)
		{
			xi = points[i].x;
			yi = points[i].y;
			zi = points[i].z;

			xj = points[j].x;
			yj = points[j].y;
			zj = points[j].z;

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

bool ComplexCreator::startCauchy(vector<int> & ligandsPermutation, string flagsFile)
{
	int sumChelation = orderAllLigands();

	if (sumChelation > maxChelation)
		return false;

	vector<double> points = getPoints(sumChelation);

	setInitialPositionCauchy(points, ligandsPermutation, flagsFile); //view fredmudar

	return true;
}



vector<double> ComplexCreator::getPoints(int totalChelation)
{
    CauchyIndex ci_(totalChelation);
	vector<CoordXYZ> mol0 = ci_.getPoints();
	size_t size = totalChelation;
	vector<double> out(3 * size);
	for (size_t i = 0; i < size; i++)
	{
		out[i] = mol0[i].x;
		out[i + size] = mol0[i].y;
		out[i + 2 * size] = mol0[i].z;
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



void ComplexCreator::setInitialPositionCauchy(
	const vector<double> &points,
	vector<int> & ligandsPermutation,
	string flagsFile)
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
		/*
		//fredmudar
		// vou entrar com -- m1(0) m2(1) m3(2) m4(3) m4(4) m4(5)
		// tenho 6 pontos -- p1 p2 p3 p4 p5 p6
		vector<int> zeroPerm(6);
		zeroPerm[0] = 1;//m1
		zeroPerm[1] = 4;//m2
		zeroPerm[2] = 3;//m3
		zeroPerm[3] = 0;//m4 
		zeroPerm[4] = 2;//b1 --- esses numeros vem do bidentateChosen
		zeroPerm[5] = 5;//b2
		// na funcao de gerar os caras eu tenho que gerar o LumpacViewInput-Permutation
		// nessa ordem m1 m2 m3 de ligantes
		*/

		// primeiro vem os assimetricos e depois os simetricos, pelo numero dos tipos eu sei.
		CauchyIndex ci_(nPoints);
		vector<int> initialPermut = ci_.zeroPermutation(flagsFile);

		int k = 0;
		for (size_t i = 0; i < allLigands.size(); i++)
		{
			vector<int> pointsOverLigand(allLigands[i].getChelation());

			// adicionar aqui uma permutacao zero para os bidentados.
			// para os monos tambem.
			// se o bi for 0-1
			// o mono esta em segundo mas a posicao dele e a 2
			// 0  1  0  0  2  0  0  0    -1  bidentates:  3  0  7  5  6  2    -1  
			// na permutacao 012345 o bidentado fica no 1 e 5 por exemplo.
			// o monodentado 1 fica na posicao 1, o monodentado 2 fica na posicao 4
			// CONFERIR - OLHAR OS GERADOS PELO PRINT MOLECULE E OS DAQUI
			// coordenada por coordenada

			for (size_t j = 0; j < allLigands[i].getChelation(); j++)
			{
				vector<int>::iterator it = find(ligandsPermutation.begin(), ligandsPermutation.end(), initialPermut[k]);
				pointsOverLigand[j] = distance(ligandsPermutation.begin(), it);
//				cout << "point" << i << "  " << pointsOverLigand[j] << endl; //fredmudar
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

