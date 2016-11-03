#include "DoAllPermutations.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "Coordstructs.h"
#include "FindIsomers.h"
#include "BuildComplex.h"
#include "RootMeanSquareDeviation.h"

using namespace std;

DoAllPermutations::DoAllPermutations() {}

DoAllPermutations::~DoAllPermutations(){}

void DoAllPermutations::calculateAll(
	string methodOptimize, 
	string methodCosmo, 
	string projectName,
	string filePermutations,
	int nPermutations)
{
	getAllPermutations(filePermutations, nPermutations);
	for (size_t i = 0; i < permutations.size(); i++)
	{
		vector<int> permutationI = permutations[i];
		string projectNameI = projectName + "-" + permutationToString(permutationI);

		buildComplexWithAPermutation(
			permutations[i],
			methodOptimize,
			methodCosmo,
			projectNameI);

	}
}

void DoAllPermutations::analysis(string crystalFile, string isomerFile)
{
	ifstream readEnergy_(isomerFile.c_str());
	string auxline, dummy;
	double energy;
	stringstream line;
	getline(readEnergy_, auxline);
	getline(readEnergy_, auxline);
	line << auxline;
	line >> dummy >> energy;
	readEnergy_.close();

	RootMeanSquareDeviation rmsd_;

	vector<CoordXYZ> crystal = rmsd_.readCoord(crystalFile);
	vector<CoordXYZ> isomer = rmsd_.readCoord(isomerFile);

	double rmsd = rmsd_.inertiaTensorComparisson(crystal, isomer);

	ofstream csv_;
	csv_.open("allEnergies", std::ofstream::out | std::ofstream::app);
	csv_ << isomerFile << " ; " << setprecision(16) << energy << " ; " << setprecision(16) << rmsd << endl;
	csv_.close();

}

void DoAllPermutations::buildComplexWithAPermutation(
	std::vector<int> permutation,
	std::string methodOptimize,
	std::string methodCosmo,
	std::string projectName
	)
{
	FindIsomers fd_;
	vector<string> options;
	string mopacExecPath;
	vector<CoordXYZ> allAtoms = fd_.buildComplexWithSelectedIsomer(
		permutation,
		projectName,
		methodOptimize,
		options,
		mopacExecPath);
	BuildComplex bc_;
	options[1] = projectName + "-water";
	options[2] = methodCosmo;
	options[4] = "";

/* 3->parametros externos 4->metal central */

	bc_.runMopacAndPrint(options, mopacExecPath, allAtoms);
}

void DoAllPermutations::getAllPermutations(string filePermutations, int nPermutations)
{
	ifstream allPermFile_(filePermutations.c_str());

	string auxline;
	getline(allPermFile_, auxline);
	while (auxline != "")
	{
		vector<int> permutationI(nPermutations);
		stringstream line1, line2;
		line1 << auxline;
		int nAtoms;
		line1 >> nAtoms;
		getline(allPermFile_, auxline);
		line2 << auxline;
		for (int i = 0; i < nPermutations; i++)
			line2 >> permutationI[i];

		permutations.push_back(permutationI);

		for (int i = 0; i < nAtoms; i++)
			getline(allPermFile_, auxline);

		getline(allPermFile_, auxline);
	}
}

string DoAllPermutations::permutationToString(vector<int>& permutation)
{
	stringstream line;
	for (size_t i = 0; i < permutation.size(); i++)
		line << permutation[i] << "-";

	string permToString;
	getline(line, permToString);
	return permToString;
}




