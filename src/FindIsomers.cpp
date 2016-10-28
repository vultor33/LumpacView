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
#include "AuxMath.h"

using namespace std;

FindIsomers::FindIsomers()
{
	identicalStructuresLimit = 0.1;
	useFile = false;
	wall0 = clock();
}

FindIsomers::~FindIsomers(){}

void FindIsomers::start()
{
	readInput();

//	fileAllIsomers = "Lumpac-allIsomers.xyz";
	if (exists_test0(fileAllIsomers))
		remove(fileAllIsomers.c_str());


	BuildComplex bc_;
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA(vector<int>(),inputInformations);
	int permutationsNumber = bc_.getLigandsPermutation().size();
	vector<CoordXYZ> atomsOriginal = ligandToCoordXYZ(allAtomsOriginal);
	vector<int> firstPermutation(permutationsNumber);
	for (int i = 0; i < permutationsNumber; i++)
		firstPermutation[i] = i;
	aplyPermutationBidentate(firstPermutation, atomsOriginal);
	if (useFile)
	{
		streamAllIsomers_.open(fileAllIsomers, std::ofstream::out | std::ofstream::app);
		appendPrintCoordXYZ(atomsOriginal, streamAllIsomers_, permutationToString(firstPermutation));
		streamAllIsomers_.close();
	}
	else
	{
		allConfigurations.push_back(atomsOriginal);
		allConfigurationPermutation.push_back(firstPermutation);
	}
	counter = 1;
	allCounter = 1;
	maxCounter = factorial(permutationsNumber);
	permutation(permutationsNumber);

	if (!useFile)
	{
		streamAllIsomers_.open(fileAllIsomers, std::ofstream::out | std::ofstream::app);
		for (size_t i = 0; i < allConfigurations.size(); i++)
			appendPrintCoordXYZ(allConfigurations[i], streamAllIsomers_, permutationToString(allConfigurationPermutation[i]));

		streamAllIsomers_ << endl << endl << "number of configurations = " << counter << endl;
		streamAllIsomers_.close();
	}
}

void FindIsomers::printSelectedIsomer(
	vector<int> permutation, 
	vector<std::string> inputInformations,
	string outputName)
{
//	if (exists_test0(outputName))
//		remove(outputName.c_str());

	appendPrintCoordXYZ(generateSelectedIsomer(permutation, inputInformations), outputName, permutationToString(permutation));
}

vector<CoordXYZ> FindIsomers::generateSelectedIsomer(vector<int> permutation, vector<string > inputInformations)
{
	BuildComplex bc_;
	vector<Ligand> allLigands = bc_.assembleComplexWithoutSA(permutation, inputInformations);
	return ligandToCoordXYZ(allLigands);
}

void FindIsomers::buildComplexWithSelectedIsomer()
{
	ifstream input_("LumpacView-inputPermutation.txt");
	int nLigands;
	string auxline;
	stringstream line1;
	getline(input_, auxline);
	line1 << auxline;
	line1 >> nLigands;
	vector<string> inputInformations(nLigands + 2);
	stringstream line2;
	getline(input_, auxline);
	line2 << auxline;
	line2 >> inputInformations[0];
	stringstream line3;
	getline(input_, auxline);
	line3 << auxline;
	line3 >> inputInformations[1];
	for (size_t i = 0; i < nLigands; i++)
	{
		stringstream linei;
		getline(input_, auxline);
		linei << auxline;
		linei >> inputInformations[i + 2];
	}
	int nPermutations;
	stringstream line4;
	getline(input_, auxline);
	line4 << auxline;
	line4 >> nPermutations;
	vector<int> permutation(nPermutations);
	for (size_t i = 0; i < nPermutations; i++)
	{
		stringstream linei;
		getline(input_, auxline);
		linei << auxline;
		linei >> permutation[i];
	}
	string execpah;
	stringstream line5;
	getline(input_, auxline);
	line5 << auxline;
	line5 >> execpah;
	string projectName;
	stringstream line6;
	getline(input_, auxline);
	line6 << auxline;
	line6 >> projectName;

	vector<string> options(5);
	options[0] = "mopac2009";
	options[1] = projectName;
	options[2] = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	options[3] = inputInformations[1];
	options[4] = inputInformations[0];

	FindIsomers fd_;
	remove("permutation-teste.xyz");
	vector<CoordXYZ> allAtoms = fd_.generateSelectedIsomer(permutation, inputInformations);
	BuildComplex bc_;
	bc_.runMopacAndPrint(options, execpah, allAtoms);



	/* EXEMPLO
	inputInformations[0] = "Eu";
	inputInformations[1] = "Eu_spk";
	inputInformations[2] = "Lumpac-View-Ligand-BUVXAR11";
	inputInformations[3] = "Lumpac-View-Ligand-BUVXAR11";
	inputInformations[4] = "Lumpac-View-Ligand-DUCNAQ-OONO";
	vector<int> permutation(4);
	permutation[0] = 0;
	permutation[1] = 1;
	permutation[2] = 2;
	permutation[3] = 3;
	vector<string> options(5);
	options[0] = "mopac2009";
	options[1] = "teste-mopac";
	options[2] = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	options[3] = inputInformations[1];
	options[4] = inputInformations[0];
	string execpah = "M2009_Ln_Orbitals.exe";
	*/
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
	allCounter++;
	RootMeanSquareDeviation rmsd_;
	BuildComplex bc_;
	vector<Ligand> allLigands = bc_.assembleComplexWithoutSA(permutation,inputInformations);
	vector<CoordXYZ> atomsPointPermutation = ligandToCoordXYZ(allLigands);

	aplyPermutationBidentate(permutation, atomsPointPermutation);
	if (atomsPointPermutation.size() == 0)
		return;
	int nMax = allLigands.size();
	/////////////////////////////////

	bool isDifferent = doOverlayWithPreviousConfigurations(atomsPointPermutation);
	if (isDifferent)
	{
		if (useFile)
		{
			streamAllIsomers_.open(fileAllIsomers, std::ofstream::out | std::ofstream::app);
			appendPrintCoordXYZ(atomsPointPermutation, streamAllIsomers_, permutationToString(permutation));
			streamAllIsomers_.close();
		}
		else
		{
			allConfigurations.push_back(atomsPointPermutation);
			allConfigurationPermutation.push_back(permutation);
		}
		counter++;
	}
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

void FindIsomers::aplyPermutationBidentate(vector<int> permutation, vector<CoordXYZ>& atomsPointPermutation)
{
	if (bidentateAtoms.size() == 0)
		return;

	int iBi1, iBi2;
	CoordXYZ meanI;
	double angle;
	AuxMath auxMath_;
	for (size_t i = 0; i < bidentateAtoms.size(); i+=2)
	{
		iBi1 = bidentateAtoms[i];
		iBi2 = bidentateAtoms[i + 1];
		angle = auxMath_.angleFrom3Points(
			atomsPointPermutation[iBi1].x, atomsPointPermutation[iBi1].y, atomsPointPermutation[iBi1].z,
			0.0e0, 0.0e0, 0.0e0,
			atomsPointPermutation[iBi2].x, atomsPointPermutation[iBi2].y, atomsPointPermutation[iBi2].z);

		if (angle > bidentateAngleCut)
		{
			vector<CoordXYZ> dummy;
			atomsPointPermutation = dummy;
			return;
		}

		meanI.x = 0.5e0*(atomsPointPermutation[iBi1].x + atomsPointPermutation[iBi2].x);
		meanI.y = 0.5e0*(atomsPointPermutation[iBi1].y + atomsPointPermutation[iBi2].y);
		meanI.z = 0.5e0*(atomsPointPermutation[iBi1].z + atomsPointPermutation[iBi2].z);
		meanI.atomlabel = "He";
		atomsPointPermutation.push_back(meanI);
	}
}


unsigned int FindIsomers::factorial(unsigned int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}

std::vector<CoordXYZ> FindIsomers::ligandToCoordXYZ(std::vector<Ligand> & allLigands)
{
	
	vector<CoordXYZ> newAllAtoms;
	/*
	(1);
	newAllAtoms[0].atomlabel = "Eu";
	newAllAtoms[0].x = 0.0e0;
	newAllAtoms[0].y = 0.0e0;
	newAllAtoms[0].z = 0.0e0;
	*/
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

void FindIsomers::appendPrintCoordXYZ(vector<CoordXYZ> & allAtoms, ofstream & fName_, string title)
{
	fName_ << allAtoms.size() << endl << title << endl;
	for (size_t i = 0; i < allAtoms.size(); i++)
	{
		fName_ << allAtoms[i].atomlabel << "  "
			<< allAtoms[i].x << "  "
			<< allAtoms[i].y << "  "
			<< allAtoms[i].z << endl;
	}
}

void FindIsomers::readInput()
{
	ifstream input_("LumpacViewInput.txt");
	string atomLabel, atomParams, auxline;
	stringstream line0, line1, line2, line3, line4, line5;
	int nAtoms, nBidentates;
	getline(input_, auxline);
	line0 << auxline;
	line0 >> fileAllIsomers;
	getline(input_, auxline);
	line1 << auxline;
	line1 >> atomLabel;
	getline(input_, auxline);
	line2 << auxline;
	line2 >> atomParams;
	getline(input_, auxline);
	line3 << auxline;
	line3 >> nAtoms;
	inputInformations.resize(nAtoms + 2);
	inputInformations[0] = atomLabel;
	inputInformations[1] = atomParams;
	for (int i = 0; i < nAtoms; i++)
	{
		stringstream linei;
		string auxRead;
		getline(input_, auxline);
		linei << auxline;
		linei >> auxRead;
		inputInformations[i + 2] = auxRead;
	}
	getline(input_, auxline);
	line4 << auxline;
	line4 >> bidentateAngleCut;
	getline(input_, auxline);
	line5 << auxline;
	line5 >> nBidentates;
	if (nBidentates > 0)
	{
		bidentateAtoms.resize(nBidentates);
		bidentateAtomsCombination.resize(nBidentates);
		for (int i = 0; i < nBidentates; i++)
		{
			stringstream linej;
			int auxReadj;
			getline(input_, auxline);
			linej << auxline;
			linej >> auxReadj;
			bidentateAtoms[i] = auxReadj;
		}
	}
	input_.close();

	AuxMath auxMath_;
	bidentateAngleCut *= (auxMath_._pi / 180);

	/*
	inputInformations[0] = "Eu";
	inputInformations[1] = "Eu_spk";
	inputInformations[2] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate1";
	inputInformations[3] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate2";
	inputInformations[4] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate3";
	inputInformations[5] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate4";
	inputInformations[6] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate5";
	inputInformations[7] = "auxLigands/Lumpac-View-Dummy-Ligand-Monodentate6";
	bidentateAngleCut = 100;
	bidentateAtoms.resize(4);
	bidentateAtomsCombination.resize(4);
	bidentateAtoms[0] = 2;
	bidentateAtoms[1] = 3;
	bidentateAtoms[2] = 4;
	bidentateAtoms[3] = 5;
	*/
}
