#include "FindIsomers.h"

#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>

#include "BuildComplex.h"
#include "Coordstructs.h"
#include "RootMeanSquareDeviation.h"

using namespace std;

FindIsomers::~FindIsomers()
{
}

FindIsomers::FindIsomers(){}

void FindIsomers::start()
{
	fileAllIsomers = "Lumpac-allIsomers.txt";
	if (exists_test0(fileAllIsomers))
		remove(fileAllIsomers.c_str());

	// H, B and C are different atoms
	// mas tambem tenho q alimentalo com a permutacao zero - mas pode ser vazio.
	// e so alimentar o assembleComplex com: vector<string> 
	// que nao precisa do input.
	// rms overlay - se o flag dos atomos for diferentes eu tenho q retornar
	// com um numero bem grande.
	// eu tenho q aplicar a permutacao
	// minimizar a combinacao dela com os atomos
	// e guardar.
	BuildComplex bc_;
	vector< string > inputInformations(6);
	inputInformations[0] = "Eu";
	inputInformations[1] = "Eu_spk";
	inputInformations[2] = "Lumpac-View-Dummy-Ligand-Monodentate";
	inputInformations[3] = "Lumpac-View-Dummy-Ligand-Monodentate";
	inputInformations[4] = "Lumpac-View-Dummy-Ligand-Monodentate";
	inputInformations[5] = "Lumpac-View-Dummy-Ligand-Monodentate";

	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA(vector<int>(),inputInformations);// permutation
	int permutationsNumber = bc_.getLigandsPermutation().size();

	streamAllIsomers_.open(fileAllIsomers, std::ofstream::out | std::ofstream::app);

	appendPrintCoordXYZ(allAtomsOriginal, fileAllIsomers, "initial configuration");
	streamAllIsomers_.close();
}



void FindIsomers::permutation(int nMax)
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

		// aplicar o vector<int> permutation

	} while (std::next_permutation(myints, myints + nMax));

	delete[] myints;
}

void FindIsomers::ligandFilePositionPermutation(vector<int> & permutation)
{
	RootMeanSquareDeviation rmsd_;
	//vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());

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
		/*
		double rmsI = rmsd_.rmsOverlay(molCrystal, molTempI);

		if (rmsI < bestRmsd)
		{
			bestRmsd = rmsI;
			bestLigand = localFilePermutation;
		}
		*/
	} while (std::next_permutation(myints, myints + nMax));

	delete[] myints;
}


vector< Ligand > FindIsomers::setThisPermutationLig(vector<int> permutation, vector<Ligand> & ligOriginal)
{
	vector<Ligand> LigPermutation(permutation.size());
	for (size_t k = 0; k < permutation.size(); k++)
		LigPermutation[k] = ligOriginal[permutation[k]];
	return LigPermutation;
}

unsigned int FindIsomers::factorial(unsigned int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}

std::vector<CoordXYZ> FindIsomers::ligandToCoordXYZ(std::vector<Ligand> & allLigands)
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

void FindIsomers::appendPrintCoordXYZ(vector<Ligand> & allAtoms, string fName, string title)
{
	appendPrintCoordXYZ(ligandToCoordXYZ(allAtoms), fName, title);
}

void FindIsomers::appendPrintCoordXYZ(vector<CoordXYZ> & allAtoms, string fName, string title)
{
	ofstream pr_;
	pr_.open (fName.c_str(), std::ofstream::out | std::ofstream::app);
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

