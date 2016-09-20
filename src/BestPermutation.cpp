#include "BestPermutation.h"

#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>

#include "Coordstructs.h"
#include "Ligand.h"
#include "RootMeanSquareDeviation.h"
#include "BuildComplex.h"

using namespace std;

BestPermutation::BestPermutation(std::string referenceFile_in)
{
	referenceFile = referenceFile_in;
	bestRmsd = 1.0e99;
}

BestPermutation::~BestPermutation(){}

void BestPermutation::findBestPermutation()
{
	BuildComplex bc_;
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA();
	int permutationsNumber = bc_.getLigandsPermutation().size();
	ligandPointPositionPermutaion(permutationsNumber);

	// PRINTING FINAL RESULT
	vector<CoordXYZ> allAtoms = ligandToCoordXYZ(bestLigand);
	stringstream convert;
	convert << bestRmsd;
	string bestRmsdTitle;
	convert >> bestRmsdTitle;
	bestRmsdTitle = "final rmsd:  " + bestRmsdTitle;
	printCoordXYZ(allAtoms, "finalConformation", bestRmsdTitle);
	////////////////////////	
}

void BestPermutation::ligandPointPositionPermutaion(int nMax)
{
	int * myints;
	myints = new int[nMax];
	for (int i = 0; i < nMax; i++)
		myints[i] = i;
	std::sort(myints, myints + nMax);
	long int size = factorial(nMax);
	vector<int> permutation(nMax);
	do
	{
		for (int i = 0; i < nMax; i++)
			permutation[i] = myints[i];

		ligandFilePositionPermutation(permutation);

	} while (std::next_permutation(myints, myints + nMax));

	delete[] myints;
}

void BestPermutation::ligandFilePositionPermutation(vector<int> & permutation)
{
	RootMeanSquareDeviation rmsd_;
	vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());

	// aplicar permutation ao sistema
	BuildComplex bc_;
	vector<Ligand> allLigands = bc_.assembleComplexWithoutSA(permutation);
	int nMax = allLigands.size();
	/////////////////////////////////

	int * myints;
	myints = new int[nMax];
	for (int i = 0; i < nMax; i++)
		myints[i] = i;
	std::sort(myints, myints + nMax);
	long int size = factorial(nMax);
	vector<int> internalPermutation(nMax);
	do
	{
		for (int i = 0; i < nMax; i++)
			internalPermutation[i] = myints[i];

		vector< Ligand > localFilePermutation = setThisPermutationLig(internalPermutation, allLigands);

		vector<CoordXYZ> molTempI = ligandToCoordXYZ(localFilePermutation);

		double rmsI = rmsd_.rmsOverlay(molCrystal, molTempI);

		if (rmsI < bestRmsd)
		{
			bestRmsd = rmsI;
			bestLigand = localFilePermutation;
		}

	} while (std::next_permutation(myints, myints + nMax));

	delete[] myints;
}

vector< Ligand > BestPermutation::setThisPermutationLig(vector<int> permutation, vector<Ligand> & ligOriginal)
{
	vector<Ligand> LigPermutation(permutation.size());
	for (size_t k = 0; k < permutation.size(); k++)
		LigPermutation[k] = ligOriginal[permutation[k]];
	return LigPermutation;
}


void BestPermutation::printCoordXYZ(vector<CoordXYZ> & allAtoms, string fName, string title)
{
	ofstream pr_(fName.c_str());
	pr_ << allAtoms.size() << endl << title << endl;
	for (size_t i = 0; i < allAtoms.size(); i++)
	{
		pr_ << allAtoms[i].atomlabel << "  "
			<< allAtoms[i].x << "  "
			<< allAtoms[i].y << "  "
			<< allAtoms[i].z << endl;
	}
	pr_.close();
}

vector<CoordXYZ> BestPermutation::ligandToCoordXYZ(vector<Ligand> & allLigands)
{
	vector<CoordXYZ> newAllAtoms(1);
	newAllAtoms[0].atomlabel = "Eu";
	newAllAtoms[0].x = 0.0e0;
	newAllAtoms[0].y = 0.0e0;
	newAllAtoms[0].z = 0.0e0;
	for (size_t i = 0; i < allLigands.size(); i++)
	{
		vector<CoordXYZ> ligandAdd = allLigands[i].getAllAtoms();
		newAllAtoms.insert(
			newAllAtoms.end(),
			ligandAdd.begin(),
			ligandAdd.end());
	}
	return newAllAtoms;
}

vector< vector<int> > BestPermutation::allFactorialPermutations(const int nMax)
{
	int * myints;
	myints = new int[nMax];
	for (int i = 0; i < nMax; i++)
		myints[i] = i;

	std::sort(myints, myints + nMax);

	long int size = factorial(nMax);

	vector< vector<int> > permutations(size);
	for (size_t i = 0; i < permutations.size(); i++)
		permutations[i].resize(nMax);

	int k = 0;
	do
	{
		permutations[k].resize(nMax);
		for (int i = 0; i < nMax; i++)
		{
			permutations[k][i] = myints[i];

			// bla bla bla

			//cout << myints[i] << "  ";
		}
		//cout << endl;
		k++;
	} while (std::next_permutation(myints, myints + nMax));

	delete[] myints;

	return permutations;
}

unsigned int BestPermutation::factorial(unsigned int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}


