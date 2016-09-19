#include "BestPermutation.h"

#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "Coordstructs.h"
#include "Ligand.h"
#include "RootMeanSquareDeviation.h"
#include "BuildComplex.h"

using namespace std;

BestPermutation::BestPermutation(vector<string> originalPermutation_in, std::string referenceFile_in)
{
	originalPermutation = originalPermutation_in;
	referenceFile = referenceFile_in;
}

BestPermutation::~BestPermutation(){}

void BestPermutation::findBestPermutation()
{
	// eficiencia eu posso ler o arquivo so um vez e guardar as informacoes na memoria.
	/*
	SAO DOIS PROBLEMAS - UM E A ORDEM DOS LIGANTES NO COMPLEXO
	ESSA DEVE SER A MESMA DO CRISTALOGRAFICO

	OUTRO E A FORMA QUE O COMPLEXO FOI MONTADO ORIGINALMENTE
	ISSO E RELEVANTE APENAS QUANDO OS LIGANTES SAO DIFERENTES
	*/

	BuildComplex bc_;

	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA();

	int permutationsNumber = bc_.getLigandsPermutation().size();

	vector< vector<int> > allPerm = allFactorialPermutations(permutationsNumber);

	vector< vector<int> > internalPerm = allFactorialPermutations(bc_.getLigandsNumber());

	RootMeanSquareDeviation rmsd_;
	vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());
	double lowestRmsd = 1.0e99;
	int lowestFilePosition = -1;
	int lowestPermutationPosition = -1;
	ofstream printPermutations_("permutations-" + referenceFile + ".txt");
	for (size_t i = 0; i < allPerm.size(); i++)
	{
		double bestRms;
		int rmsI;
		findMapToReferencePermutation(i, allPerm, internalPerm, rmsI, bestRms);
		if (bestRms < lowestRmsd)
		{
			lowestRmsd = bestRms;
			lowestFilePosition = i;
			lowestPermutationPosition = rmsI;
		}
		printPermutations_ << "iteration:" << i << " rms:" << bestRms << " --> ";
		for (size_t j = 0; j < allPerm[i].size(); j++)
			printPermutations_ << allPerm[i][j] << "  ";
		printPermutations_ << endl;
	}
	printPermutations_ << endl << endl << "best=>  " << lowestFilePosition << " : " << lowestRmsd << endl;
	printPermutations_.close();

	printSupersition(lowestFilePosition, lowestPermutationPosition, allPerm, internalPerm);
}


void BestPermutation::printSupersition(
	int lowestFilePosition, 
	int lowestPermutationPosition, 
	vector< vector<int> > & allPerm, 
	vector< vector<int> > & internalPerm)
{
	RootMeanSquareDeviation rmsd_;

	vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());

	vector<CoordXYZ> superpositionAll = molCrystal;

	BuildComplex bc_;

	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA(allPerm[lowestFilePosition]);

	vector<Ligand> allAtoms = setThisPermutationLig(internalPerm[lowestPermutationPosition], allAtomsOriginal);

	vector<CoordXYZ> molTempI = ligandToCoordXYZ(allAtoms);

	double rmsI = rmsd_.rmsOverlay(molCrystal, molTempI);

	superpositionAll.insert(
		superpositionAll.end(),
		molTempI.begin(),
		molTempI.end());

	printCoordXYZ(superpositionAll, "znormal-depois-super.xyz");
}

void BestPermutation::findMapToReferencePermutation(
	int filePermutation, 
	vector< vector<int> > & allPerm, 
	vector< vector<int> > & internalPerm, 
	int & mapToReferenceI, 
	double & mapToReferenceRms)
{
	//encontrar a configuracao que mapeia a origem
	BuildComplex bcInternal_;
	vector<Ligand> allAtomsOriginal = bcInternal_.assembleComplexWithoutSA(allPerm[filePermutation]);
	RootMeanSquareDeviation rmsd_;
	vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());
	double lowestInternalRmsd = 1.0e99;
	int lowestInernalPosition = -1;



	/* PRINTING GEOMETRIC COMBINATIONS*/
	vector<CoordXYZ> molCombinations = ligandToCoordXYZ(allAtomsOriginal);
	ofstream pr_;
	pr_.open("point-permutations.xyz", ofstream::app);
	pr_ << molCombinations.size() << endl << "useless" << endl;
	for (size_t i = 0; i < molCombinations.size(); i++)
	{
		pr_ << molCombinations[i].atomlabel << "  "
			<< molCombinations[i].x << "  "
			<< molCombinations[i].y << "  "
			<< molCombinations[i].z << endl;
	}
	pr_.close();
	/* END OF PRINTING*/



	//ofstream printPermutations_("permutations-internal" + referenceFile + ".txt");
	for (size_t i = 0; i < internalPerm.size(); i++)
	{
		vector<Ligand> allAtoms = setThisPermutationLig(internalPerm[i], allAtomsOriginal);

		vector<CoordXYZ> molTempI = ligandToCoordXYZ(allAtoms);

		double rmsI = rmsd_.rmsOverlay(molCrystal, molTempI);

		if (rmsI < lowestInternalRmsd)
		{
			lowestInernalPosition = i;
			lowestInternalRmsd = rmsI;
		}
	}
	mapToReferenceI = lowestInernalPosition;
	mapToReferenceRms = lowestInternalRmsd;

/*                           APAGAR SEM MEDO
	RootMeanSquareDeviation rmsd_;
	vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());
	double lowestRmsd = 1.0e99;
	int lowestPosition = -1;
	//ofstream printPermutations_("permutations-internal" + referenceFile + ".txt");
	for (size_t i = 0; i < allPerm.size(); i++)
	{
		BuildComplex bc_;

		vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA(allPerm[filePermutation]);

		vector<Ligand> allAtoms = setThisPermutationLig(allPerm[i], allAtomsOriginal);

		vector<CoordXYZ> molTempI = ligandToCoordXYZ(allAtoms);

		double rmsI = rmsd_.rmsOverlay(molCrystal, molTempI);

		if (rmsI < lowestRmsd)
		{
			lowestPosition = i;
			lowestRmsd = rmsI;
		}

		printPermutations_ << "iteration:" << i << " rms:" << rmsI << " --> ";
		for (size_t j = 0; j < allPerm[i].size(); j++)
			printPermutations_ << allPerm[i][j] << "  ";
		printPermutations_ << endl;
	}
	mapToReferenceI = lowestPosition;
	mapToReferenceRms = lowestRmsd;
	//printPermutations_.close();
	*/
}


vector< string > BestPermutation::setThisPermutation(vector<int> permutation)
{
	vector<string> filePermutation(originalPermutation.size());
	for (size_t k = 0; k < permutation.size(); k++)
		filePermutation[k] = originalPermutation[permutation[k]];
	return filePermutation;
}

vector< Ligand > BestPermutation::setThisPermutationLig(vector<int> permutation, vector<Ligand> & ligOriginal)
{
	vector<Ligand> LigPermutation(originalPermutation.size());
	for (size_t k = 0; k < permutation.size(); k++)
		LigPermutation[k] = ligOriginal[permutation[k]];
	return LigPermutation;
}


void BestPermutation::printCoordXYZ(vector<CoordXYZ> & allAtoms, string fName)
{
	ofstream pr_(fName.c_str());
	pr_ << allAtoms.size() << endl << "useless" << endl;
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



unsigned int BestPermutation::factorial(unsigned int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
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
