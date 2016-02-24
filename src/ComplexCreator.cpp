#include "ComplexCreator.h"

#include <vector>
#include <iostream>

#include "Ligand.h"
#include "Points.h"

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



	return true;
}

// posso seguir na ideia do Simas.
// posso estabelecer a situacao pela construcao dos pontos.

vector<double> ComplexCreator::getPoints(int totalChelation)
{
	Points allPoints;
	vector<double> out;
	switch (totalChelation)
	{
	case 2:
		out = arrayToVector(allPoints.p2, 2);
		break;
	case 3:
		out = arrayToVector(allPoints.p3, 3);
		break;
	case 4:
		out = arrayToVector(allPoints.p4, 4);
		break;
	case 5:
		out = arrayToVector(allPoints.p5, 5);
		break;
	case 6:
		out = arrayToVector(allPoints.p6, 6);
		break;
	case 7:
		out = arrayToVector(allPoints.p7, 7);
		break;
	case 8:
		out = arrayToVector(allPoints.p8, 8);
		break;
	case 9:
		out = arrayToVector(allPoints.p9, 9);
		break;
	case 10:
		out = arrayToVector(allPoints.p10, 10);
		break;
	default:
		cout << "Nao tenho essa quantidade de pontos" << endl;
		exit(1);
	}
	return out;
}

vector<double> ComplexCreator::arrayToVector(
	const double * array_in, 
	size_t size)
{
	vector<double> out(size);
	for (size_t i = 0; i < size; i++)
		out[i] = array_in[i];

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









