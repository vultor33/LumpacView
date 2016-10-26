#include "FindIsomers.h"

#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <time.h>

#include "BuildComplex.h"
#include "Coordstructs.h"
#include "RootMeanSquareDeviation.h"

using namespace std;

FindIsomers::FindIsomers()
{
	identicalStructuresLimit = 0.1;
	useFile = true;
	wall0 = clock();
}

FindIsomers::~FindIsomers(){}

void FindIsomers::start()
{
	fileAllIsomers = "Lumpac-allIsomers.xyz";
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

	//consertar o bidentado no 5
	//colocar a permutacao no titulo
	//criar uma funcao que remonte so com a permutacao.

	BuildComplex bc_;
	inputInformations.resize(10);
	inputInformations[0] = "Eu";
	inputInformations[1] = "Eu_spk";
	inputInformations[2] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate1";
	inputInformations[3] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate2";
	inputInformations[4] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate3";
	inputInformations[5] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate3";
	inputInformations[6] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate4";
	inputInformations[7] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate4";
	inputInformations[8] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate5";
	inputInformations[9] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate5";
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA(vector<int>(),inputInformations);
	int permutationsNumber = bc_.getLigandsPermutation().size();
	if (useFile)
	{
		streamAllIsomers_.open(fileAllIsomers, std::ofstream::out | std::ofstream::app);
		appendPrintCoordXYZ(allAtomsOriginal, fileAllIsomers, "initial configuration");
		streamAllIsomers_.close();
	}
	else
	{
		vector<int> firstPermutation(permutationsNumber);
		for (int i = 0; i < permutationsNumber; i++)
			firstPermutation[i] = i;
		allConfigurations.push_back(ligandToCoordXYZ(allAtomsOriginal));
		allConfigurationPermutation.push_back(firstPermutation);
	}

	counter = 1;
	allCounter = 1;
	maxCounter = factorial(permutationsNumber);
	permutation(permutationsNumber);

	cout << "number of configurations  = " << counter << endl;
}

void FindIsomers::printSelectedIsomer(
	vector<int> permutation, 
	vector<std::string> inputInformations,
	string outputName)
{
//	if (exists_test0(outputName))
//		remove(outputName.c_str());

	BuildComplex bc_;
	vector<Ligand> allLigands = bc_.assembleComplexWithoutSA(permutation, inputInformations);
	vector<CoordXYZ> atomsPointPermutation = ligandToCoordXYZ(allLigands);
	appendPrintCoordXYZ(atomsPointPermutation, outputName, permutationToString(permutation));
}

void FindIsomers::permutation(int nMax)
{
	int * myints;
	myints = new int[nMax];
	for (int i = 0; i < nMax; i++)
		myints[i] = i;
	std::sort(myints, myints + nMax);
	long int size = factorial(nMax);
	vector<int> permutationV(nMax);
	do
	{
		for (int i = 0; i < nMax; i++)
			permutationV[i] = myints[i];

		ligandFilePositionPermutation(permutationV);

	} while (std::next_permutation(myints, myints + nMax));

	delete[] myints;
}

void FindIsomers::ligandFilePositionPermutation(vector<int> & permutation)
{
	RootMeanSquareDeviation rmsd_;
	//vector<CoordXYZ> molCrystal = rmsd_.readCoord(referenceFile.c_str());
	// aplicar permutation ao sistema
	BuildComplex bc_;
	vector<Ligand> allLigands = bc_.assembleComplexWithoutSA(permutation,inputInformations);
	vector<CoordXYZ> atomsPointPermutation = ligandToCoordXYZ(allLigands);
	int nMax = allLigands.size();
	/////////////////////////////////

	bool isDifferent = doOverlayWithPreviousConfigurations(atomsPointPermutation);
	if (isDifferent)
	{
		if (useFile)
		{
			streamAllIsomers_.open(fileAllIsomers, std::ofstream::out | std::ofstream::app);
			appendPrintCoordXYZ(atomsPointPermutation, fileAllIsomers, permutationToString(permutation));
			streamAllIsomers_.close();
		}
		else
		{
			allConfigurations.push_back(atomsPointPermutation);
			allConfigurationPermutation.push_back(permutation);
		}
		counter++;
	}
	allCounter++;
	if (allCounter % 100 == 0)
	{
		wall1 = clock();
		cout << "Progress:  " << (int)(((float)allCounter / (float)maxCounter) * 100)
			<< "   elapsed time:  " << wall1 - wall0 << endl;
		wall0 = wall1;
	}

}

bool FindIsomers::doOverlayWithPreviousConfigurations(vector<CoordXYZ> & atomsPointPermutation)
{
	if (useFile)
	{
		ifstream streamAllIsomersRead_;
		streamAllIsomersRead_.open(fileAllIsomers);
		RootMeanSquareDeviation rmsd_;
		vector<CoordXYZ> atomsConfigurationsOnFile;
		//LOOP ON FILE
		do
		{
			atomsConfigurationsOnFile = readMidXyz(streamAllIsomersRead_);
			if (atomsConfigurationsOnFile.size() == 0)
				break;

			double rmsd = rmsd_.hardRmsOverlay(atomsPointPermutation, atomsConfigurationsOnFile);
			if (rmsd < identicalStructuresLimit)
				return false;

		} while (atomsConfigurationsOnFile.size() != 0);
	}
	else
	{
		RootMeanSquareDeviation rmsd_;
		vector<CoordXYZ> atomsConfigurationsOnFile;
		//LOOP ON FILE
		size_t k = 0;
		do
		{
			if (k == allConfigurations.size())
				break;

			atomsConfigurationsOnFile = allConfigurations[k];
			double rmsd = rmsd_.hardRmsOverlay(atomsPointPermutation, atomsConfigurationsOnFile);
			if (rmsd < identicalStructuresLimit)
				return false;
			else
				k++;

		} while (atomsConfigurationsOnFile.size() != 0);
	}

	return true;
}

string FindIsomers::permutationToString(vector<int>& permutation)
{
	stringstream line;
	for (size_t i = 0; i < permutation.size(); i++)
		line << permutation[i] << " ";

	string permToString;
	getline(line, permToString);
	return permToString;
}

vector<CoordXYZ> FindIsomers::readMidXyz(ifstream & openStream_)
{
	string auxline;
	getline(openStream_, auxline);
	if (auxline == "")
		return vector<CoordXYZ>();

	int natoms;
	stringstream line1;
	line1 << auxline;
	line1 >> natoms;
	vector<CoordXYZ> atomsRead(natoms);
	getline(openStream_, auxline);
	for (int i = 0; i < natoms; i++)
	{
		stringstream linei;
		getline(openStream_, auxline);
		linei << auxline;
		linei >> atomsRead[i].atomlabel
			>> atomsRead[i].x
			>> atomsRead[i].y
			>> atomsRead[i].z;
	}
	if (atomsRead.size() == 0)
	{
		cout << "error on FindIsomers::readMidXyz" << endl;
		exit(1);
	}
	else
		return atomsRead;
}

vector< Ligand > FindIsomers::setThisPermutationLig(vector<int> permutation, vector<Ligand> & ligOriginal)
{
	vector<Ligand> LigPermutation(permutation.size());
	for (size_t k = 0; k < permutation.size(); k++)
		LigPermutation[k] = ligOriginal[permutation[k]];
	return LigPermutation;
}

std::vector< CoordXYZ > FindIsomers::setThisPermutationAtoms(std::vector<int> permutation, std::vector<CoordXYZ> &  originAtoms)
{
	vector<CoordXYZ> permutAtoms(permutation.size());
	for (size_t k = 0; k < permutation.size(); k++)
		permutAtoms[k] = originAtoms[permutation[k]];
	return permutAtoms;
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



