#include "CauchyIndex.h"

#define UNIX

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <time.h>
#include <cmath>
#include <iomanip>
#ifdef UNIX
	#include <unistd.h>
#endif

#include "AllMolecularFormulas.h"
#include "RootMeanSquareDeviation.h"
#include "Coordstructs.h"
#include "AuxMath.h"

using namespace std;

CauchyIndex::CauchyIndex(int iSystem)
{
	bidentateLabels.resize(6);
	bidentateLabels[0] = "He";
	bidentateLabels[1] = "Ne";
	bidentateLabels[2] = "Ar";
	bidentateLabels[3] = "Kr";
	bidentateLabels[4] = "Xe";
	bidentateLabels[5] = "Rn";
	atomLabels.resize(12);
	atomLabels[0] = "Li";
	atomLabels[1] = "B";
	atomLabels[2] = "C";
	atomLabels[3] = "N";
	atomLabels[4] = "O";
	atomLabels[5] = "F";
	atomLabels[6] = "Cl";
	atomLabels[7] = "Br";
	atomLabels[8] = "I";
	atomLabels[9] = "Ti";
	atomLabels[10] = "Sc";
	atomLabels[11] = "Au";


	calculateAllIndexes(iSystem);
}

CauchyIndex::~CauchyIndex(){}

/*
No caso dos bidentados eu continuo mudando apenas os atomos ligantes
o negocio e que o bidentado precisa de uma funcao que determina
qual ponto preto precisa ser pintado dependendo da escolha dos atomos.
caso nenhum possa ser pintado, significa que aquele bidentado
esta com um angulo de quelação maior do que devia e também e descartado.

eu tenho que colocar os dois atomos na posicao inicial do vetor
e fazer as rotacoes de novo. no espaço dos bidentados terao espacos vazios n tem problema.

bidentados - esqueca o ponto preto, concentre-se nos atomos coordenados.

*/

std::vector<CoordXYZ> CauchyIndex::getPoints()
{
	return mol0;
}


void CauchyIndex::generateAllIndependentIsomers()
{
	systemSize = mol0.size();
	long int size = factorial(systemSize); // 10+ explosion

	nRotations = allRotationTransforms.size() + 1;
	vector<liPermutation> allPermutations;

	int * myints;
	myints = new int[systemSize];
	for (size_t i = 0; i < systemSize; i++)
		myints[i] = i;
	std::sort(myints, myints + systemSize);
	vector<int> permutation(systemSize);
	bool equal;
	int k = 0; // fredmudar
	do
	{
		k++;
		if (k % 10000 == 0)
			cout << "k:  " << k << endl;

		equal = false;
		for (size_t i = 0; i < systemSize; i++)
			permutation[i] = myints[i];

		// check if this permutation was already studied
		for (size_t i = 0; i < allPermutations.size(); i++)
		{
			for (size_t j = 0; j < nRotations; j++)
			{
				if (permutation == allPermutations[i].rotPermutations[j])
				{
					equal = true;
					break;
				}
				if (equal)
					break;
			}
		}

		if (equal)
			continue;
		else
		{
			//add permutation
			vector< vector<int> > rotPermutations(nRotations);
			rotPermutations[0] = permutation;
			for (size_t j = 1; j < nRotations; j++)
				rotPermutations[j] = applyRotation(permutation, j - 1);
			liPermutation auxLiPermutation;
			auxLiPermutation.rotPermutations = rotPermutations;
			allPermutations.push_back(auxLiPermutation);
		}
	} while (std::next_permutation(myints, myints + systemSize));

	ofstream of_("printFinalPermutations.txt");

	for (size_t i = 0; i < allPermutations.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << allPermutations[i].rotPermutations[0][j] << "  ";
		}
		of_ << endl;
	}

	delete[] myints;
}

void CauchyIndex::generateAllIndependentIsomersIO()
{
	systemSize = mol0.size();
	long int size = factorial(systemSize); // 10+ explosion
	nRotations = allRotationTransforms.size() + 1;

	int * myints;
	myints = new int[systemSize];
	for (size_t i = 0; i < systemSize; i++)
		myints[i] = i;
	std::sort(myints, myints + systemSize);
	vector<int> permutation(systemSize);
	remove("allCauchyRotatations.txt");
	string cauchyFileName = "a llCauchyRotatations.txt";
	bool equal;
	int k = 0;

	clock_t begin = clock();
	do
	{
		equal = false;
		for (size_t i = 0; i < systemSize; i++)
			permutation[i] = myints[i];

		// check if this permutation was already studied
		if (k != 0)
		{
			ifstream allPermutOpen_(cauchyFileName.c_str());
			while(true)
			{
				vector<int> permutationK = readCauchyNotations(allPermutOpen_);
				if (permutationK.size() == 0)
					break;
				if (permutation == permutationK)
				{
					equal = true;
					break;
				}
			}
		}

		if (equal)
			continue;
		else
		{
			//add permutation
			vector< vector<int> > rotPermutations(nRotations);
			rotPermutations[0] = permutation;
			for (size_t j = 1; j < nRotations; j++)
				rotPermutations[j] = applyRotation(permutation, j - 1);
			writeCauchyRotations(cauchyFileName, rotPermutations);
		}

		if (k % 10000 == 0)
			cout << "k:  " << (double)k / (double)size << endl;
		k++;
	} while (std::next_permutation(myints, myints + systemSize));

	delete[] myints;

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "demorou:  " << elapsed_secs << "  segundos" << endl;
}


void CauchyIndex::generateAllIndependentIsomersRuntimeRotations()
{
	systemSize = mol0.size();
	int size = factorial(systemSize); // 10+ explosion
	nRotations = allRotationTransforms.size() + 1;
	vector< vector<int> > allPermutations(size / nRotations);
	size_t isomerCounter = 0;
	bool equal;

	int * myints;
	myints = new int[systemSize];
	for (size_t i = 0; i < systemSize; i++)
		myints[i] = i;
	std::sort(myints, myints + systemSize);
	vector<int> permutation(systemSize);

	int k = 0; //fredmudar
	do
	{
		k++;
		if (k % 1000 == 0)
			cout << "k:  " << k << endl;

		equal = false;

		for (size_t i = 0; i < systemSize; i++)
			permutation[i] = myints[i];

		for (size_t i = 0; i < isomerCounter; i++)
		{
			equal = compareTwoIsomers(allPermutations[i], permutation);
			if (equal)
				break;
		}
		if (!equal)
		{
			allPermutations[isomerCounter] = permutation;
			isomerCounter++;
		}

	} while (std::next_permutation(myints, myints + systemSize));

	ofstream of_("printFinalPermutations.txt");

	for (size_t i = 0; i < allPermutations.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << allPermutations[i][j] << "  ";
		}
		of_ << endl;
	}

	delete[] myints;
}


void CauchyIndex::generateAllIndependentIsomersRuntimeRotationsAndReadBlock(string blockFileName)
{
	systemSize = mol0.size();
	int size = factorial(systemSize);
	nRotations = allRotationTransforms.size() + 1;
	vector< vector<int> > allPermutations;
	vector< vector<int> > blockPermutations;
	ifstream openedFile_(blockFileName.c_str());
	bool equal;
	while (true)
	{
		vector<int> auxBlock = readCauchyNotations(openedFile_);
		if (auxBlock.size() == 0)
			break;

		equal = false;

		for (size_t i = 0; i < allPermutations.size(); i++)
		{
			equal = compareTwoIsomers(allPermutations[i], auxBlock);
			if (equal)
				break;
		}
		if (!equal)
			allPermutations.push_back(auxBlock);
	}

/*  SAME CODE WITH RAM
	while (true)
	{
		vector<int> auxBlock = readCauchyNotations(openedFile_);
		if (auxBlock.size() == 0)
			break;
		blockPermutations.push_back(auxBlock);
	}
	for (size_t iBlock = 0; iBlock < blockPermutations.size(); iBlock++)
	{
		equal = false;

		for (size_t i = 0; i < allPermutations.size(); i++)
		{
			equal = compareTwoIsomersAtoms(allPermutations[i], blockPermutations[iBlock]);
			if (equal)
				break;
		}
		if (!equal)
			allPermutations.push_back(blockPermutations[iBlock]);
	}
*/

	ofstream of_(("isomersOf-" + blockFileName).c_str());
	for (size_t i = 0; i < allPermutations.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << allPermutations[i][j] << "  ";
		}
		of_ << endl;
	}
}


void CauchyIndex::correctIndependentIsomers()
{
	int systemSize = 12;
	ifstream openedFile_("independent-isomers-1");
	ofstream out_("independent");
	while (true)
	{
		vector<int> auxBlock2 = readCauchyNotations(openedFile_);
		if (auxBlock2.size() == 0)
			break;
		vector<int> auxBlock(systemSize);
		auxBlock[0] = 0;
		for (size_t i = 0; i < systemSize - 1; i++)
			auxBlock[i + 1] = auxBlock2[i] + 1;

		for (size_t i = 0; i < auxBlock.size(); i++)
			out_ << auxBlock[i] << " ";
		out_ << endl;
	}	
	openedFile_.close();
	out_.close();
}


void CauchyIndex::generateAllIndependentIsomers12(string blockFileName)
{
        systemSize = mol0.size();
        int size = factorial(systemSize);
        nRotations = allRotationTransforms.size() + 1;
        vector< vector<int> > allPermutations;
        vector< vector<int> > blockPermutations;
        ifstream openedFile_(blockFileName.c_str());
        bool equal;
        while (true)
        {
                vector<int> auxBlock2 = readCauchyNotations(openedFile_);
                if (auxBlock2.size() == 0)
                        break;
		vector<int> auxBlock(systemSize);
		auxBlock[0] = 0;
		for(size_t i = 0; i < systemSize - 1; i++)
			auxBlock[i+1] = auxBlock2[i] + 1;

                equal = false;

                for (size_t i = 0; i < allPermutations.size(); i++)
                {
                        equal = compareTwoIsomers(allPermutations[i], auxBlock);
                        if (equal)
                                break;
                }
                if (!equal)
                        allPermutations.push_back(auxBlock);
        }
        ofstream of_(("isomersOf-" + blockFileName).c_str());
        for (size_t i = 0; i < allPermutations.size(); i++)
        {
                for (size_t j = 0; j < systemSize; j++)
                {
                        of_ << allPermutations[i][j] << "  ";
                }
                of_ << endl;
        }
}


void CauchyIndex::generateAllIndependentIsomersWithFlag(
	string blockFileName, 
	string flagsFile, 
	string code)
{
	vector<int> atomTypes;
	vector<int> bidentateAtomChosen;
	readAtomTypesAndBidentateChosenFile(flagsFile, atomTypes, bidentateAtomChosen);

//	molecularFormulaToCauchyCode(code, atomTypes, bidentateAtomChosen);

	systemSize = mol0.size();
	int size = factorial(systemSize);
	nRotations = allRotationTransforms.size() + 1;
	vector< vector<int> > allPermutations;
	vector< vector<int> > blockPermutations;
	ifstream openedFile_(blockFileName.c_str());
	int compareResult;
	bool equal;
	while (true)
	{
		vector<int> auxBlock = readCauchyNotations(openedFile_);
		if (auxBlock.size() == 0)
			break;

		equal = false;

		for (size_t i = 0; i < allPermutations.size(); i++)
		{
			compareResult = compareTwoIsomersWithLabelsRotations(
				atomTypes,
				bidentateAtomChosen,
				allPermutations[i],
				auxBlock);

			if (compareResult == 0 || compareResult == 2)
				equal = true;

			if (equal)
				break;
		}
		if (!equal)
			allPermutations.push_back(auxBlock);
	}

	ofstream of_((code + "-" + blockFileName).c_str());
/* HEADER
	of_ << allPermutations.size() << "  types:  ";
	for (size_t i = 0; i < atomTypes.size(); i++)
		of_ << atomTypes[i] << "  ";
	of_ << "  -1  bidentates:  ";
	for (size_t i = 0; i < bidentateAtomChosen.size(); i++)
		of_ << bidentateAtomChosen[i] << "  ";
	of_ << "  -1  ";
	of_ << endl;
*/
	for (size_t i = 0; i < allPermutations.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << allPermutations[i][j] << "  ";
		}
		of_ << endl;
	}

	/*print xyz
	ofstream xyzFile_((code + "-" + blockFileName + ".xyz").c_str());
	for (size_t i = 0; i < allPermutations.size(); i++)
	{
	printMolecule(
	allPermutations[i],
	atomTypes,
	bidentateAtomChosen,
	xyzFile_);
	}
	xyzFile_.close();
	*/


	/* TESTE ESPECIFICO
	vector<int> auxBlock(6);
	auxBlock[0] = 0;
	auxBlock[1] = 1;
	auxBlock[2] = 4;
	auxBlock[3] = 2;
	auxBlock[4] = 3;
	auxBlock[5] = 5;
	equal = false;
	int i = 2;
	compareResult = compareTwoIsomersWithLabels(
	atomTypes,
	bidentateAtomChosen,
	allPermutations[i],
	auxBlock);
	*/
}


vector<int> CauchyIndex::zeroPermutation(string flagsFile)
{
	vector<int> atomTypes;
	vector<int> bidentateAtomChosen;
	readAtomTypesAndBidentateChosenFile(flagsFile, atomTypes, bidentateAtomChosen);

	//os ligantes vao estar na ordem -- C, B e m.

	vector<int> zero(atomTypes.size());
	int k = 0;
	for (size_t i = 0; i < bidentateAtomChosen.size(); i++)
	{
		zero[k] = bidentateAtomChosen[i];
		k++;
	}
	for (size_t i = 0; i < atomTypes.size(); i++)
	{
		if (find(bidentateAtomChosen.begin(), bidentateAtomChosen.end(), i) != bidentateAtomChosen.end())
			atomTypes[i] = 999;
	}
	while (k != (int)atomTypes.size())
	{
		vector<int>::iterator lowest = min_element(atomTypes.begin(), atomTypes.end());
		int lowestIndex = distance(atomTypes.begin(), lowest);
		zero[k] = lowestIndex;
		k++;
		atomTypes[lowestIndex] = 999;
	}
	return zero;
}




void CauchyIndex::molecularFormulaToCauchyCode(
	string code,
	vector<int> & atomTypes,
	vector<int> & bidentateAtomsChosen)
{
	AllMolecularFormulas allMolecular_;
	atomTypes.resize(mol0.size());
	vector< vector<int> > molecularCode = allMolecular_.stringToNumber(code);
	int kType = 0;

	int bidentateProblem;
	vector<int> permutation(mol0.size());
	for (size_t i = 0; i < mol0.size(); i++)
		permutation[i] = i;
    int k = 0;
	do
	{
	    k++;
	    if(k > 100000)
        {
            cout << "CauchyIndex::molecularFormulaToCauchyCode error" << endl;
            exit(1);
        }
		vector<int> auxBidentateAtomsChosen;
		vector<int> auxAtomTypes = atomTypes;
		int kAuxBidentatePosition = 0;
		int kTypes = 0;
		for (size_t i = 0; i < molecularCode[2].size(); i++)
		{
			for (size_t j = 0; j < molecularCode[2][i]; j++)
			{
				setBidentateChosen(auxBidentateAtomsChosen);
				auxAtomTypes[auxBidentateAtomsChosen[kAuxBidentatePosition]] = kTypes;
				kAuxBidentatePosition++;
				auxAtomTypes[auxBidentateAtomsChosen[kAuxBidentatePosition]] = kTypes + 1;
				kAuxBidentatePosition++;
			}
			kTypes += 2;
		}
		for (size_t i = 0; i < molecularCode[1].size(); i++)
		{
			for (size_t j = 0; j < molecularCode[1][i]; j++)
			{
				setBidentateChosen(auxBidentateAtomsChosen);
				auxAtomTypes[auxBidentateAtomsChosen[kAuxBidentatePosition]] = kTypes;
				kAuxBidentatePosition++;
				auxAtomTypes[auxBidentateAtomsChosen[kAuxBidentatePosition]] = kTypes;
				kAuxBidentatePosition++;
			}
            kTypes++;
		}
		vector<int> auxAtomsChosen = auxBidentateAtomsChosen;
        for (size_t i = 0; i < molecularCode[0].size(); i++)
        {
            for (size_t j = 0; j < molecularCode[0][i]; j++)
            {
                int iAuxAtomTypes;
                do
                {
                    iAuxAtomTypes = auxMath_.randomNumber(0, mol0.size() - 1);
                } while (find(auxAtomsChosen.begin(), auxAtomsChosen.end(), iAuxAtomTypes) != auxAtomsChosen.end());
                auxAtomTypes[iAuxAtomTypes] = kTypes;
                auxAtomsChosen.push_back(iAuxAtomTypes);
            }
            kTypes++;
        }
		bidentateProblem = compareTwoIsomersWithLabelsRotations(
			atomTypes,
			auxBidentateAtomsChosen,
			permutation,
			permutation);

		if (bidentateProblem != 2)
		{
			bidentateAtomsChosen = auxBidentateAtomsChosen;
			atomTypes = auxAtomTypes;
			break;
		}
	} while (true);

	//fredapagar
	/*
	vector<string> labels(6);
	labels[0] = "F";
	labels[1] = "B";
	labels[2] = "C";
	labels[3] = "Li";
	labels[4] = "N";
	labels[5] = "O";
	vector<string> atomLabels(6);
	//tenho q ter uma coisa que pega o types e transforma em letras
	for(size_t i = 0; i < atomTypes.size(); i++)
    {
            atomLabels[i] = labels[atomTypes[i]];
            cout << "labels:  "  << atomTypes[i] << endl;
    }
    for(size_t i = 0; i < bidentateAtomsChosen.size(); i++)
        cout << "bidentate:  " << bidentateAtomsChosen[i] << endl;

    ofstream of_("teste-print.xyz");
    printMolecule(
		permutation,
		atomLabels,
		bidentateAtomsChosen,
		of_);
    of_.close();
	*/
	//for (size_t i = 0; i < mol0.size(); i++)
	//	cout << permutation[i] = i;
}


void CauchyIndex::setBidentateChosen(std::vector<int>& bidentateAtomsChosen)
{
	int i1, i2;
	do
	{
		i1 = auxMath_.randomNumber(0, mol0.size() - 1);
	} while (find(bidentateAtomsChosen.begin(), bidentateAtomsChosen.end(), i1) != bidentateAtomsChosen.end());
	bidentateAtomsChosen.push_back(i1);
	do
	{
		i2 = auxMath_.randomNumber(0, mol0.size() - 1);
	} while (find(bidentateAtomsChosen.begin(), bidentateAtomsChosen.end(), i2) != bidentateAtomsChosen.end());
	bidentateAtomsChosen.push_back(i2);
}


int CauchyIndex::compareTwoIsomersWithLabelsRotations(
	std::vector<int>& atomTypes,
	std::vector<int>& bidentateAtomsChosen,
	std::vector<int>& permutationIsomer1,
	std::vector<int>& permutationIsomer2)
{
	size_t size = mol0.size();
	vector<int> types1(size);
	vector<int> types2(size);
	for (size_t i = 0; i < size; i++)
	{
		types1[i] = atomTypes[permutationIsomer1[i]];
		types2[i] = atomTypes[permutationIsomer2[i]];
	}
	vector<int> bidPos1 = applyPermutationBidentates(permutationIsomer1, bidentateAtomsChosen);
	if (bidPos1.size() == 0 && bidentateAtomsChosen.size() != 0)
	{
		return 2;
	}
	vector<int> bidPos2 = applyPermutationBidentates(permutationIsomer2, bidentateAtomsChosen);
	if (bidPos2.size() == 0 && bidentateAtomsChosen.size() != 0)
	{
		return 2;
	}
	sort(bidPos1.begin(), bidPos1.end());
	sort(bidPos2.begin(), bidPos2.end());
	if (types1 == types2)
	{
		if (bidPos1 == bidPos2)
			return 0;
	}
 	for (size_t j = 0; j < allRotationTransforms.size(); j++)
	{
		vector<int> auxPerm = applyRotation(permutationIsomer2, j);
		for (size_t i = 0; i < size; i++)
		{
			//types1[i] = atomTypes[permutationIsomer1[i]]; fredapagar
			types2[i] = atomTypes[auxPerm[i]];
		}
		if (types1 == types2)
		{
			vector<int> bidPos2 = applyPermutationBidentates(auxPerm, bidentateAtomsChosen);
			sort(bidPos2.begin(), bidPos2.end());
			if (bidPos2.size() == 0 && bidentateAtomsChosen.size() != 0)
			{
				//fred apagar -- rotacoes nunca levam bidentados a posicoes proibidas.
				cout << "CauchyIndex::compareTwoIsomersWithLabels( stop" << endl;
				exit(1);
			}
			if (bidPos1 == bidPos2)
				return 0;
		}
	}
	return 1;
}


int CauchyIndex::compareTwoIsomersWithLabels(
	std::vector<int>& atomTypes,
	std::vector<int>& bidentateAtomsChosen,
	std::vector<int>& permutationIsomer1,
	std::vector<int>& permutationIsomer2)
{
	size_t size = mol0.size();
	vector<int> types1(size);
	vector<int> types2(size);
	for (size_t i = 0; i < size; i++)
	{
		types1[i] = atomTypes[permutationIsomer1[i]];
		types2[i] = atomTypes[permutationIsomer2[i]];
	}
	vector<int> bidPos1;
	vector<int> bidPos2;
	if (bidentateAtomsChosen.size() != 0)
	{
		bidPos1 = applyPermutationBidentates(permutationIsomer1, bidentateAtomsChosen);
		bidPos2 = applyPermutationBidentates(permutationIsomer2, bidentateAtomsChosen);
		if (bidPos2.size() == 0)
			return 2;
		sort(bidPos1.begin(), bidPos1.end());
		sort(bidPos2.begin(), bidPos2.end());
	}
	if (types1 == types2)
	{
		if (bidPos1 == bidPos2)
			return 0;
	}
	return 1;
}


bool CauchyIndex::compareTwoIsomers(
	std::vector<int>& permutationIsomer1,
	std::vector<int>& permutationIsomer2)
{
	size_t nRotations = allRotationTransforms.size();
	if(mol0.size() == 12)
		nRotations = 4;

	for (size_t j = 0; j < nRotations; j++)
	{
		vector<int> auxPerm = applyRotation(permutationIsomer2, j);
		if (permutationIsomer1 == auxPerm)
			return true;
	}
	return false;
}


void CauchyIndex::mergeBlocks(std::vector<std::string> & allBlockNames, int nPieces)
{
	remove("auxAllBlocks.txt");
	int total = 0;
	for(size_t i = 0; i < allBlockNames.size(); i++)
	{
		ifstream openedFile_(allBlockNames[i].c_str());
		while(true)
		{
			vector<int> auxBlock = readCauchyNotations(openedFile_);
                	if (auxBlock.size() == 0)
                        	break;
			printCauchyNotation("auxAllBlocks.txt",auxBlock);
			total++;
		}
		openedFile_.close();
	}
	cout << "total permutations:  " << total << endl;

	ifstream openedFile_("auxAllBlocks.txt");
	int cutBlock = total / nPieces;
	int iBlockNumber = 0;
	int iBlock = 0;
	string fileName = "new-block--";
	string blockFileName;
	while(true)
	{
		vector<int> auxBlock = readCauchyNotations(openedFile_);
                if (auxBlock.size() == 0)
			break;
		if (iBlock % cutBlock == 0)
		{
                        iBlockNumber++;
                        stringstream convert;
                        convert << iBlockNumber;
                        string auxBlockName;
                        convert >> auxBlockName;
                        blockFileName = fileName + auxBlockName + ".txt";
                        remove(blockFileName.c_str());
	    	}
                printCauchyNotation(blockFileName, auxBlock);
                iBlock++;
	}
	openedFile_.close();

}


void CauchyIndex::printBlock(int nPieces)
{
	int systemSize = mol0.size();
	int totalSize = factorial(systemSize);
	int cutBlock = totalSize / nPieces;
	int * myints;
	myints = new int[systemSize];
	for (size_t i = 0; i < systemSize; i++)
		myints[i] = i;
	std::sort(myints, myints + systemSize);
	vector<int> permutation(systemSize);

	stringstream convert0;
	convert0 << systemSize;
	string auxFileName;
	convert0 >> auxFileName;
	string fileName = "block-" + auxFileName + "---";
	string blockFileName;
	int iBlock = 0;
	int iBlockNumber = 0;
	do
	{
		for (size_t i = 0; i < systemSize; i++)
			permutation[i] = myints[i];

		if (iBlock % cutBlock == 0)
		{
			iBlockNumber++;
			stringstream convert;
			convert << iBlockNumber;
			string auxBlockName;
			convert >> auxBlockName;
			blockFileName = fileName + auxBlockName + ".txt";
			remove(blockFileName.c_str());
		}

		printCauchyNotation(blockFileName, permutation);
		iBlock++;
	} while (std::next_permutation(myints, myints + systemSize));

	delete[] myints;

}


void CauchyIndex::generateBlockFiles(int n, int kInit, int kFinal)
{
        int * myints;
        myints = new int[n];
        for (size_t i = 0; i < n; i++)
                myints[i] = i;
	string permutName;
	int k = 1;
	do
        {
		if(k >= kInit && k <= kFinal)
		{
			stringstream convert;
	                for(int i = 0; i < n; i++)
        	                convert << myints[i] << "-";
			permutName = convert.str();

			ofstream name_(permutName.c_str());
			name_.close();
		}
		k++;

        } while (std::next_permutation(myints, myints + n));

	delete[] myints;

}


void CauchyIndex::generateRAMBlock(int n, int kInit, int kFinal, vector< vector<int> > & ramBlock)
{
	ramBlock.resize(kFinal - kInit + 1);
        int * myints;
        myints = new int[n];
        for (size_t i = 0; i < n; i++)
                myints[i] = i;
        int k = 1;
	int ramPos = 0;
	vector<int> auxPerm(n);
        do
        {
                if(k >= kInit && k <= kFinal)
                {
			for(size_t i = 0; i < n; i++)
				auxPerm[i] = myints[i];
			ramBlock[ramPos] = auxPerm;
			ramPos++;
                }
                k++;

        } while (std::next_permutation(myints, myints + n));
        delete[] myints;
}


void CauchyIndex::generateSlurmFilesToDeletion(int nSystem, int nProc, string machineType)
{
	int deletionSystem;
	if (nSystem == 12)
	{
		deletionSystem = 12;
		nSystem--;
	}
	else
		deletionSystem = nSystem;

	int totalSize = factorial(nSystem);
	if(totalSize % (nProc * nProc) != 0)
	{
		cout << "system and proc don't agree" << endl;
		exit(1);
	}
	int bigBlock = totalSize / nProc;
	int smallBlock = bigBlock / nProc;

	ofstream fileRunAll_("runBlock.x");
	fileRunAll_ << "#!/bin/bash" << endl;
	for(int i = 1; i <= nProc; i++)
	{
		stringstream convert3;
		convert3 << i;
		if(machineType == "slurm")
			fileRunAll_ << "sbatch " << convert3.str() << ".srm" << endl;
		else if(machineType == "pc")
			fileRunAll_ << "./" << convert3.str() << ".x" << endl;
		else
		{
			cout << "machine type not found" << endl;
			exit(1);
		}
	}
	system("chmod u+x runBlock.x");

	// 2 -- nProc
	string blockName = "Block";
	for (int i = 2; i <= nProc; i++)
	{
		stringstream convert;
		convert << i;
		system(("mkdir " + blockName + convert.str()).c_str());
		system(("cp lumpacview.exe  " + blockName + convert.str()).c_str());

		int indexBlock = (i - 1)*bigBlock + 1;
		int jSmallInit, jSmallEnd;
		for (int j = 1; j <= nProc; j++)
		{
			jSmallInit = indexBlock + (j - 1)*smallBlock;
			jSmallEnd = indexBlock + j*smallBlock - 1;
			stringstream convert2;
			convert2 << j;
			if (machineType == "slurm")
			{
				ofstream fileSrm_((convert2.str() + ".srm").c_str());
				fileSrm_ << "#!/bin/bash" << endl;
				fileSrm_ << "#SBATCH -n 1" << endl;
				fileSrm_ << "#SBATCH --hint=nomultithread" << endl;
				fileSrm_ << "RODADIR=/home/guga/PERMUTACOES/running/" + blockName + convert.str() << endl;
				fileSrm_ << "cd $RODADIR" << endl;
				fileSrm_ << "./lumpacview.exe  wholeBlockDeletion   " << deletionSystem << "  1  " << indexBlock - 1
					<< "  " << jSmallInit << "   " << jSmallEnd << endl;
				fileSrm_.close();
				system(("mv  " + convert2.str() + ".srm" + "  " + blockName + convert.str()).c_str());
				system(("cp runBlock.x  " + blockName + convert.str()).c_str());
			}
			else if (machineType == "pc")
			{
				ofstream fileSrm_((convert2.str() + ".x").c_str());
				fileSrm_ << "#!/bin/bash" << endl;
				fileSrm_ << "./lumpacview.exe  wholeBlockDeletion  " << deletionSystem << "  1  " << indexBlock - 1
					<< "  " << jSmallInit << "   " << jSmallEnd << endl;
				fileSrm_.close();
				system(("chmod u+x  " + convert2.str() + ".x").c_str());
				system(("mv  " + convert2.str() + ".x" + "  " + blockName + convert.str()).c_str());
				system(("cp runBlock.x  " + blockName + convert.str()).c_str());

			}
		}
	}

	// generate block 1 here.
	vector< vector<int> > ramBlock;
        generateRAMBlock(
                nSystem,
                1,
                bigBlock,
                ramBlock);
	string blockFileName1 = "independent-isomers-1";
	for(size_t i = 0; i < ramBlock.size(); i++)
	{
		printCauchyNotation(blockFileName1, ramBlock[i]);
	}
}

void CauchyIndex::generateSlurmFilesToDeletionFlags(
	int deletionSystem,
	int total,
	int bigBlockSize,
	int smallBlockSize,
	string compositionFile,
	string rawIsomersFile,
	string workingDir,
	string machineType)
{
	int nWholeBlocks = floor((double)total / (double)bigBlockSize);
	int lastBlock = bigBlockSize * nWholeBlocks + 1;
	int smallBlockNumber = bigBlockSize / smallBlockSize;

	ofstream fileRunAll_("runBlock.x");
	fileRunAll_ << "#!/bin/bash" << endl;
        for(int i = 0; i < smallBlockNumber; i++)
        {
                stringstream convert3;
                convert3 << i;
                if(machineType == "slurm")
                        fileRunAll_ << "sbatch " << convert3.str() << ".srm" << endl;
                else if(machineType == "pc")
                        fileRunAll_ << "./" << convert3.str() << ".x" << endl;
                else
                {
                        cout << "machine type not found" << endl;
                        exit(1);
                }
        }
        system("chmod u+x runBlock.x");
	fileRunAll_.close();




	// mid blocks
	string blockName = "Block";
        for (int i = nWholeBlocks; i > 1; i--)
        {
                stringstream convert;
                convert << i;
                system(("mkdir " + blockName + convert.str()).c_str());
                system(("cp lumpacview.exe  " + blockName + convert.str()).c_str());
		int indexBlock = (i-1) * bigBlockSize; 
                int jSmallInit, jSmallEnd;
                for (int j = 0; j < smallBlockNumber; j++)
                {
			jSmallInit = (i - 1) * bigBlockSize + j * smallBlockSize + 1;
			jSmallEnd = (i - 1) * bigBlockSize + (j + 1) * smallBlockSize;
                        stringstream convert2;
                        convert2 << j;
                        if (machineType == "slurm")
                        {
                                ofstream fileSrm_((convert2.str() + ".srm").c_str());
                                fileSrm_ << "#!/bin/bash" << endl;
                                fileSrm_ << "#SBATCH -n 1" << endl;
                                fileSrm_ << "#SBATCH --hint=nomultithread" << endl;
                                fileSrm_ << "RODADIR=" + workingDir + "/" + blockName + convert.str() << endl;
                                fileSrm_ << "cd $RODADIR" << endl;
                                fileSrm_ << "./lumpacview.exe  compositionBlockDeletion   " 
					<< deletionSystem << "  "
					<< rawIsomersFile << "  "
					<< compositionFile << "  "
					<< "  1  " << indexBlock
                                        << "  " << jSmallInit << "   " << jSmallEnd << endl;
                                fileSrm_.close();
                                system(("mv  " + convert2.str() + ".srm" + "  " + blockName + convert.str()).c_str());
                                system(("cp runBlock.x  " + blockName + convert.str()).c_str());
                        }
                        else if (machineType == "pc")
                        {
                                ofstream fileSrm_((convert2.str() + ".x").c_str());
                                fileSrm_ << "#!/bin/bash" << endl;
                                fileSrm_ << "./lumpacview.exe  compositionBlockDeletion   "
                                        << deletionSystem << "  "
                                        << rawIsomersFile << "  "
                                        << compositionFile << "  "
                                        << "  1  " << indexBlock
                                        << "  " << jSmallInit << "   " << jSmallEnd << endl;
                                fileSrm_.close();
                                system(("chmod u+x  " + convert2.str() + ".x").c_str());
                                system(("mv  " + convert2.str() + ".x" + "  " + blockName + convert.str()).c_str());
                                system(("cp runBlock.x  " + blockName + convert.str()).c_str());
			}
		}

	}


	int i = nWholeBlocks;
	stringstream convert;
        convert << i + 1;
        system(("mkdir " + blockName + convert.str()).c_str());
        system(("cp lumpacview.exe  " + blockName + convert.str()).c_str());
        int indexBlock = i * bigBlockSize;
        int jSmallInit, jSmallEnd;
	bool haveToBreak = false;
        for (int j = 0; j < smallBlockNumber; j++)
        {
		jSmallInit = lastBlock + j * smallBlockSize;
        	if (total > lastBlock + (j + 1) * smallBlockSize - 1)
	        {
                        jSmallEnd = lastBlock + (j + 1) * smallBlockSize - 1;
		}
		else
		{
                        jSmallEnd = total;
			haveToBreak = true;
		}

                stringstream convert2;
                convert2 << j;
                if (machineType == "slurm")
                {
                	ofstream fileSrm_((convert2.str() + ".srm").c_str());
                	fileSrm_ << "#!/bin/bash" << endl;
                        fileSrm_ << "#SBATCH -n 1" << endl;
                        fileSrm_ << "#SBATCH --hint=nomultithread" << endl;
                        fileSrm_ << "RODADIR=/home/guga/PERMUTACOES/running/" + blockName + convert.str() << endl;
                        fileSrm_ << "cd $RODADIR" << endl;
                        fileSrm_ << "./lumpacview.exe  compositionBlockDeletion   "
                                        << deletionSystem << "  "
                                        << rawIsomersFile << "  "
                                        << compositionFile << "  "
                                        << "  1  " << indexBlock
                                        << "  " << jSmallInit << "   " << jSmallEnd << endl;
                        fileSrm_.close();
                        system(("mv  " + convert2.str() + ".srm" + "  " + blockName + convert.str()).c_str());
                        system(("cp runBlock.x  " + blockName + convert.str()).c_str());
                }
                else if (machineType == "pc")
                {
                	ofstream fileSrm_((convert2.str() + ".x").c_str());
                	fileSrm_ << "#!/bin/bash" << endl;
                        fileSrm_ << "./lumpacview.exe  compositionBlockDeletion   "
                                        << deletionSystem << "  "
                                        << rawIsomersFile << "  "
                                        << compositionFile << "  "
                                        << "  1  " << indexBlock
                                        << "  " << jSmallInit << "   " << jSmallEnd << endl;
                        fileSrm_.close();
                        system(("chmod u+x  " + convert2.str() + ".x").c_str());
                        system(("mv  " + convert2.str() + ".x" + "  " + blockName + convert.str()).c_str());
                        system(("cp runBlock.x  " + blockName + convert.str()).c_str());
                }
		if(haveToBreak)
			break;
	}


	ofstream isomers1_("independent-isomers-1");
        ifstream openedFile_(rawIsomersFile.c_str());
	for(int i = 1; i <= bigBlockSize; i++)
	{
                vector<int> permutation = readCauchyNotations(openedFile_);
		for(size_t j = 0; j < permutation.size(); j++)
		{
			isomers1_ << permutation[j] << "  ";
		}
		isomers1_ << endl;
	}
	openedFile_.close();
	isomers1_.close();	

}



void CauchyIndex::runall(int blockInit, int blockFinal,string machineType,string workingDirectory)
{
#ifdef UNIX
	string folderName;
	for(int i = blockInit; i <= blockFinal; i++)
	{
                stringstream convert;
                convert << i;
                folderName = "Block" + convert.str();
                chdir(folderName.c_str());
                system("./runBlock.x");
		chdir(workingDirectory.c_str());
		cout << i << "   finished" << endl;
	}
#endif
}

void CauchyIndex::cleanBlocksAndGenerateIsomers(
	int nProc, 
	int systemSize, 
	string composition,
	string workingDirectory,
	string machineType)
{
#ifdef UNIX
	string folderName;
	for(int i = 2; i <= nProc; i++)
	{
        stringstream convert;
        convert << i;
        folderName = "Block" + convert.str();
        chdir(folderName.c_str());
        system(("cat block* > independent-isomers-" + convert.str()).c_str());
        system(("mv independent-isomers-" + convert.str() + " .. ").c_str());
		chdir(workingDirectory.c_str());
		system(("rm -rf " + folderName).c_str());
	}

	stringstream convert0;
	convert0 << systemSize;
	string systemSizeName = convert0.str();
	for(int i = 1; i <= nProc; i++)
	{
        stringstream convert;
        convert << i;

		if (machineType == "slurm")
		{
			ofstream fileSrm_((convert.str() + ".srm").c_str());
			fileSrm_ << "#!/bin/bash" << endl;
			fileSrm_ << "#SBATCH -n 1" << endl;
			fileSrm_ << "#SBATCH --hint=nomultithread" << endl;
			fileSrm_ << "RODADIR=" << workingDirectory << endl;
			fileSrm_ << "cd $RODADIR" << endl;
			fileSrm_ << ("./lumpacview.exe blockGeneration " + composition + "  " + systemSizeName + " independent-isomers-" + convert.str()).c_str()
				<< endl;
			fileSrm_.close();
			system(("sbatch " + convert.str() + ".srm").c_str());
		}
		else
			system(("./lumpacview.exe blockGeneration " + composition + "  " + systemSizeName + " independent-isomers-" + convert.str()).c_str());
	}
#endif
}



void CauchyIndex::doBlockDeletion(
	int kInit,
	int kFinal)
{
	int systemSize = mol0.size();
        int * myints;
        myints = new int[systemSize];
        for (size_t i = 0; i < systemSize; i++)
                myints[i] = i;
        string permutName;
        int k = 1;
        do
        {
                if(k >= kInit && k <= kFinal)
                {
			if (k % 100 == 0)
                        	cout << "k:  " << k << endl;
			vector<int> permutation(systemSize);
        		for (size_t i = 0; i < systemSize; i++)
        		        permutation[i] = myints[i];
			for (size_t j = 0; j < allRotationTransforms.size(); j++)
			{
				vector<int> auxPerm = applyRotation(permutation, j);
				stringstream convert;
				for(int i = 0; i < systemSize; i++)                
					convert << auxPerm[i] << "-";          	
				permutName = convert.str();
#ifdef UNIX
				unlink(permutName.c_str());
#endif
			}
                }
                k++;

        } while (std::next_permutation(myints, myints + systemSize));

        delete[] myints;
}


void CauchyIndex::doBlockRAMDeletion(
    int kInit,
    int kFinal,
	int ramInit,
	int ramFinal)
{
	int systemSize = mol0.size();
	vector< vector<int> > ramBlock;
 	generateRAMBlock(
		systemSize,
		ramInit,
		ramFinal,
		ramBlock);

	vector<int> permutation(systemSize);
        for (size_t i = 0; i < systemSize; i++)
                permutation[i] = i;
        int k = 1;
        do
        {
                if(k >= kInit && k <= kFinal)
                {
                        for (size_t j = 0; j < allRotationTransforms.size(); j++)
                        {
                                vector<int> auxPerm = applyRotation(permutation, j);
				/* OTHER DELETION METHOD
				for(size_t k = 0; k < ramBlock.size(); k++)
               			{
                        		if(auxPerm == ramBlock[k])
                        		{
                               			vector< vector<int> >::iterator it = ramBlock.begin() + k;
                                		rotate(it, it+1,ramBlock.end());
                                		ramBlock.pop_back();
                                		break;
                        		}

		                }
				*/
				ramBlock.erase(std::remove(ramBlock.begin(), ramBlock.end(), auxPerm), ramBlock.end());
                        }
                }
                k++;

        } while (std::next_permutation(permutation.begin(), permutation.end()));

	stringstream convert;
	convert << "block-" << ramInit << "-to-" << ramFinal << ".txt";
	string blockDelName = convert.str();
        ofstream of_(blockDelName.c_str());
        for (size_t i = 0; i < ramBlock.size(); i++)
        {
                for (size_t j = 0; j < systemSize; j++)
                {
                        of_ << ramBlock[i][j] << "  ";
                }
                of_ << endl;
        }
	of_.close();
}



void CauchyIndex::doBlockRAMDeletion12(
	int kInit,
	int kFinal,
	int ramInit,
	int ramFinal)
{
	int systemSize = 11;
	vector< vector<int> > ramBlock;
	generateRAMBlock(
		systemSize,
		ramInit,
		ramFinal,
		ramBlock);

	for (size_t i = 0; i < ramBlock.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			ramBlock[i][j] += 1;
		}
		ramBlock[i].insert(ramBlock[i].begin(), 0);
	}

	vector<int> permutation(systemSize);
	for (size_t i = 0; i < systemSize; i++)
		permutation[i] = i;
	int k = 1;
	do
	{
		if (k >= kInit && k <= kFinal)
		{
			vector<int> permutation12 = permutation;
			for (size_t j = 0; j < systemSize; j++)
			{
				permutation12[j] += 1;
			}
			permutation12.insert(permutation12.begin(), 0);
			for (size_t j = 0; j < 4; j++) //rotations over 0
			{
				vector<int> auxPerm = applyRotation(permutation12, j);
				ramBlock.erase(std::remove(ramBlock.begin(), ramBlock.end(), auxPerm), ramBlock.end());
			}
		}
		k++;

	} while (std::next_permutation(permutation.begin(), permutation.end()));

	stringstream convert;
	convert << "block-" << ramInit << "-to-" << ramFinal << ".txt";
	string blockDelName = convert.str();
	ofstream of_(blockDelName.c_str());
	for (size_t i = 0; i < ramBlock.size(); i++)
	{
		for (size_t j = 0; j < systemSize + 1; j++)
		{
			of_ << ramBlock[i][j] << "  ";
		}
		of_ << endl;
	}
	of_.close();
}


void CauchyIndex::generateAtomTypesAndBidentateChosenFile(string complexCode)
{
	vector<int> atomTypes;
	vector<int> bidentateChosen;
	molecularFormulaToCauchyCode(complexCode, atomTypes, bidentateChosen);
	ofstream of_((complexCode + "---atomTypes.txt").c_str());
	for (size_t i = 0; i < atomTypes.size(); i++)
		of_ << atomTypes[i] << "  ";
	of_ << "  -1  bidentates:  ";
	for (size_t i = 0; i < bidentateChosen.size(); i++)
		of_ << bidentateChosen[i] << "  ";
	of_ << "  -1  ";
	of_ << endl;
	of_.close();
}



void CauchyIndex::readAtomTypesAndBidentateChosenFile(
	std::string fileName, 
	vector<int> & atomTypes, 
	vector<int> & bidentateChosen)
{
	ifstream in_(fileName.c_str());
	string line;
	getline(in_, line);
	stringstream auxline;
	auxline << line;
	while (true)
	{
		int aux;
		auxline >> aux;
		if (aux == -1)
			break;
		else
			atomTypes.push_back(aux);
	}
	string dummy;
	auxline >> dummy;
	while (true)
	{
		int aux;
		auxline >> aux;
		if (aux == -1)
			break;
		else
			bidentateChosen.push_back(aux);
	}
	in_.close();
}


// ATENCAO  --- ativar o [[[generateAtomTypesAndBidentateChosenFile]]]
//              antes de comecar e colcoar os atomTypes no flagsFile.
void CauchyIndex::doBlockDeletionFlags(
	string skeletonFile,
	string flagsFile,
	int kInit,
	int kFinal,
	int ramInit,
	int ramFinal)
{
	int systemSize = mol0.size();
	ifstream openedFile_(skeletonFile.c_str());
	vector< vector<int> > ramBlock = readCauchyNotationsRAMBlock(openedFile_, ramInit, ramFinal);
	openedFile_.close();

	vector<int> atomTypes;
	vector<int> bidentateChosen;
	readAtomTypesAndBidentateChosenFile(flagsFile, atomTypes, bidentateChosen);

	openedFile_.open(skeletonFile.c_str());
	int k = 1;
	do
	{
		vector<int> permutation = readCauchyNotations(openedFile_);
		if (permutation.size() == 0)
			break;
		if (k >= kInit && k <= kFinal) // from block use this to delete
		{
			int compare;
			int kRam = 0;
			while (kRam < ramBlock.size())
			{

				// sem rotacoes int compare = compareTwoIsomersWithLabels
				compare = compareTwoIsomersWithLabels(
					atomTypes,
					bidentateChosen,
					permutation,
					ramBlock[kRam]);
				if ((compare == 0) || (compare == 2))
				{
					vector< vector<int> >::iterator it = ramBlock.begin() + kRam;
					rotate(it, it + 1, ramBlock.end());
					ramBlock.pop_back();
					kRam--;
				}
				kRam++;
			}
			for (size_t j = 0; j < allRotationTransforms.size(); j++)
			{
				vector<int> auxPerm = applyRotation(permutation, j);
				kRam = 0;
				while(kRam < ramBlock.size())
				{
					// sem rotacoes int compare = compareTwoIsomersWithLabels
					compare = compareTwoIsomersWithLabels(
						atomTypes,
						bidentateChosen,
						auxPerm,
						ramBlock[kRam]);
					if ((compare == 0) || (compare == 2))
					{
						vector< vector<int> >::iterator it = ramBlock.begin() + kRam;
						rotate(it, it + 1, ramBlock.end());
						ramBlock.pop_back();
						kRam--;
					}
					kRam++;
				}
			}
		}
		k++;
	} while (true);
	openedFile_.close();

	stringstream convert;
	convert << "block-" << ramInit << "-to-" << ramFinal << ".txt";
	string blockDelName = convert.str();
	ofstream of_(blockDelName.c_str());
	for (size_t i = 0; i < ramBlock.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << ramBlock[i][j] << "  ";
		}
		of_ << endl;
	}
	of_.close();
}

/*
void CauchyIndex::doBlockDeletionFlags(
string skeletonFile,
string flagsFile,
int kInit,
int kFinal,
int ramInit,
int ramFinal)
{
int systemSize = mol0.size();
ifstream openedFile_(skeletonFile.c_str());
vector< vector<int> > ramBlock = readCauchyNotationsRAMBlock(openedFile_, ramInit, ramFinal);
openedFile_.close();

vector<int> atomTypes;
vector<int> bidentateChosen;
readAtomTypesAndBidentateChosenFile(flagsFile, atomTypes, bidentateChosen);

openedFile_.open(skeletonFile.c_str());
int k = 1;
do
{
vector<int> permutation = readCauchyNotations(openedFile_);
if (permutation.size() == 0)
break;
if (k >= kInit && k <= kFinal)
{
int compare;
compare = compareTwoIsomersWithLabels(
atomTypes,
bidentateChosen,
permutation,
ramBlock[k]);
if (compare == 0)
{
vector< vector<int> >::iterator it = ramBlock.begin() + k;
rotate(it, it + 1, ramBlock.end());
ramBlock.pop_back();
}
for (size_t j = 0; j < allRotationTransforms.size(); j++)
{
vector<int> auxPerm = applyRotation(permutation, j);
for (size_t k = 0; k < ramBlock.size(); k++)
{
	// sem rotacoes int compare = compareTwoIsomersWithLabels
	compare = compareTwoIsomersWithLabels(
		atomTypes,
		bidentateChosen,
		auxPerm,
		ramBlock[k]);
	if (compare == 0)
	{
		vector< vector<int> >::iterator it = ramBlock.begin() + k;
		rotate(it, it + 1, ramBlock.end());
		ramBlock.pop_back();
	}
}
//ramBlock.erase(std::remove(ramBlock.begin(), ramBlock.end(), auxPerm), ramBlock.end()); fredapagar
			}
		}
		k++;
	} while (true);
	openedFile_.close();

	stringstream convert;
	convert << "block-" << ramInit << "-to-" << ramFinal << ".txt";
	string blockDelName = convert.str();
	ofstream of_(blockDelName.c_str());
	for (size_t i = 0; i < ramBlock.size(); i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << ramBlock[i][j] << "  ";
		}
		of_ << endl;
	}
	of_.close();
}
*/





























































void CauchyIndex::rotationTest(
	vector<string> & atoms,
	vector<int> & bidentateAtomsChosen)
{
	size_t size = mol0.size();
	vector<int> permutation(size);
	for (size_t i = 0; i < size; i++)
		permutation[i] = i;
	if (atoms.size() == 0)
	{
		atoms.resize(size);
		for (size_t i = 0; i < size; i++)
			atoms[i] = " H ";
	}
	atoms[0] = " C ";
	atoms[1] = " N ";
	atoms[2] = " Li ";
	ofstream of_("printTestRotations.xyz");
	for (size_t i = 0; i < allRotationTransforms.size(); i++)
	{
		printMolecule(permutation, atoms, bidentateAtomsChosen, of_);
		vector<int> auxPerm = applyRotation(permutation, i);
		printMolecule(auxPerm, atoms, bidentateAtomsChosen, of_);
	}
	of_.close();
}

void CauchyIndex::calculateAllIndexes(int iSystem)
{
	setSystem(iSystem);
	for (size_t i = 0; i < allRotationsMatrix.size(); i++)
	{
		vector<int> cauchyI = calculateRotationTransform(i);

		allRotationTransforms.push_back(cauchyI);
	}
	calculateBidentateMap();
}

std::vector<int> CauchyIndex::calculateRotationTransform(int rotation)
{
	vector<CoordXYZ> mol = mol0;
	vector< vector<double> > rot = allRotationsMatrix[rotation].mRot;
	size_t size = mol.size();
	vector<CoordXYZ> molRot = mol;
	for (size_t i = 0; i < size; i++)
	{
		vector<double> newCoord = auxMath_.matrixXVector(
			rot,
			mol[i].x,
			mol[i].y,
			mol[i].z);
		molRot[i].x = newCoord[0];
		molRot[i].y = newCoord[1];
		molRot[i].z = newCoord[2];
	}

// fredapagar
//	molRot[0].atomlabel = "Li";
//	mol[0].atomlabel = "Li";
//	molRot[1].atomlabel = "F";
//	mol[1].atomlabel = "F";
//	printMoleculeFast(mol);
//	printMoleculeFast(molRot);
	vector<int> cauchy(size);
	double atomPosition;
	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			atomPosition = sqrt(
				(molRot[i].x - mol[j].x) * (molRot[i].x - mol[j].x) +
				(molRot[i].y - mol[j].y) * (molRot[i].y - mol[j].y) +
				(molRot[i].z - mol[j].z) * (molRot[i].z - mol[j].z));

			//cout << "i:  " << i << " j: " << j << "  dist: " << atomPosition << endl;
			if (atomPosition < 1.0e-2)
			{
				cauchy[i] = j;
				break;
			}
			if (j == (size - 1))
			{
				cout << "rotation problem" << endl;
				exit(1);
			}
		}
	}

	//printCauchyNotation(cauchy);

	return cauchy;
}

void CauchyIndex::printCauchyNotation(vector<int> & cauchyList)
{
	size_t size = cauchyList.size();
	for (size_t i = 0; i < size; i++)
	{
		cout << i << "   ";
	}
	cout << endl;
	for (size_t i = 0; i < size; i++)
	{
		cout << cauchyList[i] << "   ";
	}
	cout << endl;
}

void CauchyIndex::printCauchyNotation(
	std::string fileName,
	std::vector<int> & cauchyList)
{
	ofstream cauchyFile_;
	cauchyFile_.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);
	for (size_t j = 0; j < cauchyList.size(); j++)
	{
		cauchyFile_ << cauchyList[j] << " ";
	}
	cauchyFile_ << endl;
	cauchyFile_.close();
}


void CauchyIndex::printMoleculeFast(std::vector<CoordXYZ>& mol)
{
	ofstream of_;
	of_.open("printFastMolecule.xyz", std::ofstream::out | std::ofstream::app);
	size_t size = mol.size();
	of_ << size << endl << "title " << endl;
	for (size_t i = 0; i < size; i++)
	{
		if (mol[i].atomlabel == "")
			mol[i].atomlabel = "H  ";

		of_ << mol[i].atomlabel << "  "
			<<  setprecision(12) << mol[i].x << "  "
			<<  setprecision(12) << mol[i].y << "  "
			<<  setprecision(12) << mol[i].z << endl;

	}
	of_.close();
}

void CauchyIndex::setAllRotations(const vector<double> & allRotationsVector)
{
	allRotationsMatrix.resize(allRotationsVector.size() / 4);
	for (size_t i = 0; i < allRotationsVector.size(); i += 4)
	{
		allRotationsMatrix[i / 4].mRot = auxMath_.rotationMatrix(
			allRotationsVector[i],
			allRotationsVector[i + 1],
			allRotationsVector[i + 2],
			allRotationsVector[i + 3]);
	}
}

void CauchyIndex::calculateBidentateMap()
{
	double xi, yi, zi, xj, yj, zj;
	int kMolBidentate = 0;
	bidentateMap.resize(mol0.size());
	for (size_t ii = 0; ii < mol0.size(); ii++)
		bidentateMap[ii].resize(mol0.size());
	for (size_t i = 0; i < mol0.size() - 1; i++)
	{
		bidentateMap[i][i] = -1;
		for (size_t j = i + 1; j < mol0.size(); j++)
		{
			xi = mol0[i].x;
			yi = mol0[i].y;
			zi = mol0[i].z;

			xj = mol0[j].x;
			yj = mol0[j].y;
			zj = mol0[j].z;

			double angle = auxMath_.angleFrom3Points(xi, yi, zi, 0.0e0, 0.0e0, 0.0e0, xj, yj, zj);

			if (angle < cutAngle)
			{
				CoordXYZ tempBidentate;
				tempBidentate.atomlabel = "C ";
				tempBidentate.x = 0.5e0 * (xi + xj);
				tempBidentate.y = 0.5e0 * (yi + yj);
				tempBidentate.z = 0.5e0 * (zi + zj);
				molBidentate.push_back(tempBidentate);
				bidentateMap[i][j] = kMolBidentate;
				bidentateMap[j][i] = kMolBidentate;
				kMolBidentate++;
			}
			else
			{
				bidentateMap[i][j] = -1;
				bidentateMap[j][i] = -1;
			}
		}
	}
}

vector<int> CauchyIndex::applyRotation(
	const vector<int> & permutation,
	int iRotation)
{
	size_t size = permutation.size();
	vector<int> permRot(size);
	for (size_t i = 0; i < size; i++)
		permRot[i] = permutation[allRotationTransforms[iRotation][i]];

	return permRot;
}

std::vector<int> CauchyIndex::applyPermutationCoordinates(
	const std::vector<int> & permutation,
	const std::vector<std::string> & atoms,
	const std::vector<int> & bidentateAtomsChosen)
{
	size_t size = atoms.size();
	for (size_t i = 0; i < size; i++)
	{
		mol0[i].atomlabel = atoms[permutation[i]];
	}
	vector<int> bidentateAtomsChosenRotated = bidentateAtomsChosen;
	for (size_t i = 0; i < bidentateAtomsChosenRotated.size(); i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			if (permutation[j] == bidentateAtomsChosenRotated[i])
			{
				bidentateAtomsChosenRotated[i] = j;
				break;
			}
		}
	}
	return bidentateAtomsChosenRotated;
}

std::vector<int> CauchyIndex::applyPermutationBidentates(
	const std::vector<int>& permutation,
	const std::vector<int>& bidentateAtomsChosen)
{
	size_t size = permutation.size();
	vector<int> bidentateAtomsChosenRotated = bidentateAtomsChosen;
	for (size_t i = 0; i < bidentateAtomsChosenRotated.size(); i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			if (permutation[j] == bidentateAtomsChosenRotated[i])
			{
				bidentateAtomsChosenRotated[i] = j;
				break;
			}
		}
	}
	vector<int> bidentatePositions;
	for (size_t i = 0; i < bidentateAtomsChosen.size(); i += 2)
	{
		int mapPositions = bidentateMap[bidentateAtomsChosenRotated[i]][bidentateAtomsChosenRotated[i + 1]];
		if (mapPositions == -1)
			return vector<int>();
		bidentatePositions.push_back(mapPositions);
	}
	return bidentatePositions;
}



void CauchyIndex::printMolecule(
	vector<int> & permutation,
	vector<int> & atomTypes,
	const vector<int> & bidentateAtomsChosen,
	ofstream & printFile_)
{
	vector<string> labels(permutation.size());
	for (size_t i = 0; i < atomTypes.size(); i++)
		labels[i] = atomLabels[atomTypes[i]];
	printMolecule(
		permutation,
		labels,
		bidentateAtomsChosen,
		printFile_
	);
}

void CauchyIndex::printMolecule(
	vector<int> & permutation,
	const vector<string> & atoms,
	const vector<int> & bidentateAtomsChosen,
	ofstream & printFile_)
{
	if (permutation.size() == 0)
	{
		permutation.resize(atoms.size());
		for (size_t i = 0; i < atoms.size(); i++)
			permutation[i] = i;
	}
	vector<int> bidentateAtomsChosenRotated = applyPermutationCoordinates(permutation, atoms, bidentateAtomsChosen);
	vector<CoordXYZ> bidentates;
	int k = 0;
	for (size_t i = 0; i < bidentateAtomsChosen.size(); i+=2)
	{
		int mapPosition = bidentateMap[bidentateAtomsChosenRotated[i]][bidentateAtomsChosenRotated[i + 1]];
		if (mapPosition == -1)
			cout << "CauchyIndex::printMolecule - wrong bidentate choice" << endl;

		CoordXYZ auxBidentateAtom = molBidentate[mapPosition];
		auxBidentateAtom.atomlabel = bidentateLabels[k];
		k++;
		bidentates.push_back(auxBidentateAtom);
	}
	int nAtoms = mol0.size() + bidentates.size();
	printFile_ << nAtoms << endl << "title" << endl;
	printFile_,setprecision(12);
	for (int i = 0; i < nAtoms; i++)
	{
		if (i < (int)mol0.size())
		{
			printFile_ << mol0[i].atomlabel << "  "
				 << fixed << setprecision(12) << setw(16) << mol0[i].x << "  "
				 << fixed << setprecision(12) << setw(16) << mol0[i].y << "  "
				<< fixed << setprecision(12) << setw(16) << mol0[i].z << endl;
		}
		else
		{
			printFile_ << bidentates[i - mol0.size()].atomlabel << "  "
				<< fixed << setprecision(12) << setw(16) << bidentates[i - mol0.size()].x << "  "
				<< fixed << setprecision(12) << setw(16) << bidentates[i - mol0.size()].y << "  "
				<< fixed << setprecision(12) << setw(16) << bidentates[i - mol0.size()].z << endl;
		}
	}
}


void CauchyIndex::printMoleculeFromFile(string fileName)
{
	ifstream allPermt_(fileName.c_str());
	string line;
	getline(allPermt_, line);
	stringstream convert;
	convert << line;
	string dummy;
	convert >> dummy;
	convert >> dummy;
	vector<int> atomTypes;
	while (true)
	{
		int type;
		convert >> type;
		if (type == -1)
			break;
		atomTypes.push_back(type);
	}
	convert >> dummy;
	vector<int> bidentateChosen;
	while (true)
	{
		int bid;
		convert >> bid;
		if (bid == -1)
			break;
		bidentateChosen.push_back(bid);
	}

	ofstream xyz_((fileName + ".xyz").c_str());
	while(true)
	{
		vector<int> permutation = readCauchyNotations(allPermt_);
		if (permutation.size() == 0)
			break;
		printMolecule(permutation, atomTypes, bidentateChosen, xyz_);
	}
	xyz_.close();
	allPermt_.close();
}





void CauchyIndex::printAllMoleculesFromFile(string composition)
{
	vector<int> atomTypes;
	vector<int> bidentateChosen;
	readAtomTypesAndBidentateChosenFile("results-" + composition, atomTypes, bidentateChosen);
	
	ifstream filePermutations_(("results-" + composition).c_str());
	string line;
	getline(filePermutations_,line);	
	while(true)
	{
		vector<int> permutation = readCauchyNotations(filePermutations_);
		if(permutation.size() == 0)
			break;
		string permtString = permutationToString(permutation);
		ofstream fileXyz_((composition + "-" + permtString + ".xyz").c_str());
		printMolecule(permutation,atomTypes,bidentateChosen,fileXyz_);
		fileXyz_.close();	
	}
	filePermutations_.close();
#ifdef UNIX
	system(("mkdir " + composition).c_str());
	system(("mv " + composition + "* " + composition).c_str());
#endif

}




unsigned int CauchyIndex::factorial(unsigned int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}

void CauchyIndex::writeCauchyRotations(string fileName, vector< vector<int> > & rotPermutations)
{
	ofstream cauchyFile_;
	cauchyFile_.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

	for (size_t i = 0; i < nRotations; i++)
	{
		for (size_t j = 0; j < systemSize; j++)
		{
			cauchyFile_ << rotPermutations[i][j] << " ";
		}
		cauchyFile_ << endl;
	}
	cauchyFile_.close();
}

vector<int> CauchyIndex::readCauchyNotations(ifstream & openendFile_)
{
	int size = mol0.size();
	vector<int> notation;
	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;
	notation.resize(size);
	stringstream line;
	line << auxline;
	for (size_t i = 0; i < size; i++)
	{
		line >> notation[i];
	}
	return notation;
}

vector< vector<int> > CauchyIndex::readCauchyNotationsRAMBlock(ifstream & openendFile_, int kInit, int kFinal)
{
	vector< vector<int> > ramBlock;
	int k = 1;
	while (true)
	{
		vector<int> isomer = readCauchyNotations(openendFile_);
		if (isomer.size() == 0)
			break;
		if (k >= kInit && k <= kFinal)
		{
			ramBlock.push_back(isomer);
		}
		k++;
	}
	return ramBlock;
}


string CauchyIndex::permutationToString(vector<int> permutation)
{
	stringstream permt;
	for(size_t i = 0; i < permutation.size(); i++)
		permt << permutation[i] << "-";	
	return permt.str();
}




















































void CauchyIndex::setSystem(int system)
{
	vector<double> vectorRotations;
	vector<double> auxReferenceAxis;

	switch (system)
	{
	case 5:
		mol0.resize(5);
		mol0[0].x = 0.000;
		mol0[0].y = 0.000;
		mol0[0].z = 1.000;
		mol0[1].x = 1.000;
		mol0[1].y = 0.000;
		mol0[1].z = 0.000;
		mol0[2].x = -1.000;
		mol0[2].y = 0.000;
		mol0[2].z = 0.000;
		mol0[3].x = 0.000;
		mol0[3].y = 0.86602540;
		mol0[3].z = -0.50000000;
		mol0[4].x = 0.000;
		mol0[4].y = -0.86602540;
		mol0[4].z = -0.50000000;
		//c3 - 1
		vectorRotations.resize(20);
		vectorRotations[0] = 1.0e0;
		vectorRotations[1] = 0.0e0;
		vectorRotations[2] = 0.0e0;
		vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
		//c3 - 2
		vectorRotations[4] = 1.0e0;
		vectorRotations[5] = 0.0e0;
		vectorRotations[6] = 0.0e0;
		vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;
		//c2 - 1
		vectorRotations[8] = mol0[0].x;
		vectorRotations[9] = mol0[0].y;
		vectorRotations[10] = mol0[0].z;
		vectorRotations[11] = auxMath_._pi;
		//c2 - 2
		vectorRotations[12] = mol0[3].x;
		vectorRotations[13] = mol0[3].y;
		vectorRotations[14] = mol0[3].z;
		vectorRotations[15] = auxMath_._pi;
		//c2 - 3
		vectorRotations[16] = mol0[4].x;
		vectorRotations[17] = mol0[4].y;
		vectorRotations[18] = mol0[4].z;
		vectorRotations[19] = auxMath_._pi;
		//cut angle
		cutAngle = 3.0e0 * auxMath_._pi / 4.0e0;
		break;


	case 6:
		mol0.resize(6);
		mol0[0].x = 0.00000000;
		mol0[0].y = 0.00000000;
		mol0[0].z = 1.00000000;
		mol0[1].x = 1.00000000;
		mol0[1].y = 0.00000000;
		mol0[1].z = 0.00000000;
		mol0[2].x = 0.00000000;
		mol0[2].y = 1.00000000;
		mol0[2].z = 0.00000000;
		mol0[3].x = -1.00000000;
		mol0[3].y = 0.00000000;
		mol0[3].z = 0.00000000;
		mol0[4].x = 0.00000000;
		mol0[4].y = -1.00000000;
		mol0[4].z = 0.00000000;
		mol0[5].x = 0.00000000;
		mol0[5].y = 0.00000000;
		mol0[5].z = -1.00000000;
		cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;
		auxReferenceAxis.resize(3);
		vectorRotations.resize(92);
		// 1 c3 por face.
		// tem 3 pontos que geram os c4 e c3.
		// 6 c2 cruzam as 12 arestas.
		//c4 - 1 (ref 0)
		vectorRotations[0] = mol0[0].x;
		vectorRotations[1] = mol0[0].y;
		vectorRotations[2] = mol0[0].z;
		vectorRotations[3] = auxMath_._pi / 2.0e0;
		//c4 - 2 (ref 0)
		vectorRotations[4] = mol0[0].x;
		vectorRotations[5] = mol0[0].y;
		vectorRotations[6] = mol0[0].z;
		vectorRotations[7] = auxMath_._pi;
		//c4 - 3 (ref 0)
		vectorRotations[8] = mol0[0].x;
		vectorRotations[9] = mol0[0].y;
		vectorRotations[10] = mol0[0].z;
		vectorRotations[11] = 3.0e0 * auxMath_._pi / 2.0e0;

		//c4 - 1 (ref 1)
		vectorRotations[12] = mol0[1].x;
		vectorRotations[13] = mol0[1].y;
		vectorRotations[14] = mol0[1].z;
		vectorRotations[15] = auxMath_._pi / 2.0e0;
		//c4 - 2 (ref 1)
		vectorRotations[16] = mol0[1].x;
		vectorRotations[17] = mol0[1].y;
		vectorRotations[18] = mol0[1].z;
		vectorRotations[19] = auxMath_._pi;
		//c4 - 3 (ref 1)
		vectorRotations[20] = mol0[1].x;
		vectorRotations[21] = mol0[1].y;
		vectorRotations[22] = mol0[1].z;
		vectorRotations[23] = 3.0e0 * auxMath_._pi / 2.0e0;

		//c4 - 1 (ref 2)
		vectorRotations[24] = mol0[2].x;
		vectorRotations[25] = mol0[2].y;
		vectorRotations[26] = mol0[2].z;
		vectorRotations[27] = auxMath_._pi / 2.0e0;
		//c4 - 2 (ref 2)
		vectorRotations[28] = mol0[2].x;
		vectorRotations[29] = mol0[2].y;
		vectorRotations[30] = mol0[2].z;
		vectorRotations[31] = auxMath_._pi;
		//c4 - 3 (ref 2)
		vectorRotations[32] = mol0[2].x;
		vectorRotations[33] = mol0[2].y;
		vectorRotations[34] = mol0[2].z;
		vectorRotations[35] = 3.0e0 * auxMath_._pi / 2.0e0;

		//c2 - 0 - 1
		auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[36] = auxReferenceAxis[0];
		vectorRotations[37] = auxReferenceAxis[1];
		vectorRotations[38] = auxReferenceAxis[2];
		vectorRotations[39] = auxMath_._pi;

		//c2 - 0 - 2
		auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[2].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[2].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[40] = auxReferenceAxis[0];
		vectorRotations[41] = auxReferenceAxis[1];
		vectorRotations[42] = auxReferenceAxis[2];
		vectorRotations[43] = auxMath_._pi;

		//c2 - 0 - 3
		auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[3].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[3].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[44] = auxReferenceAxis[0];
		vectorRotations[45] = auxReferenceAxis[1];
		vectorRotations[46] = auxReferenceAxis[2];
		vectorRotations[47] = auxMath_._pi;

		//c2 - 0 - 4
		auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[4].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[4].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[48] = auxReferenceAxis[0];
		vectorRotations[49] = auxReferenceAxis[1];
		vectorRotations[50] = auxReferenceAxis[2];
		vectorRotations[51] = auxMath_._pi;

		//c2 - 1 - 2
		auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[52] = auxReferenceAxis[0];
		vectorRotations[53] = auxReferenceAxis[1];
		vectorRotations[54] = auxReferenceAxis[2];
		vectorRotations[55] = auxMath_._pi;

		//c2 - 1 - 4
		auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[4].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[4].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[56] = auxReferenceAxis[0];
		vectorRotations[57] = auxReferenceAxis[1];
		vectorRotations[58] = auxReferenceAxis[2];
		vectorRotations[59] = auxMath_._pi;

		//c3 ---- 0 - 1 - 2
		auxReferenceAxis[0] = (mol0[0].x + mol0[1].x + mol0[2].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[1].y + mol0[2].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[1].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[60] = auxReferenceAxis[0];
		vectorRotations[61] = auxReferenceAxis[1];
		vectorRotations[62] = auxReferenceAxis[2];
		vectorRotations[63] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[64] = auxReferenceAxis[0];
		vectorRotations[65] = auxReferenceAxis[1];
		vectorRotations[66] = auxReferenceAxis[2];
		vectorRotations[67] = 4.0e0 * auxMath_._pi / 3.0e0;

		//c3 ---- 0 - 2 - 3
		auxReferenceAxis[0] = (mol0[0].x + mol0[2].x + mol0[3].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[2].y + mol0[3].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[2].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[68] = auxReferenceAxis[0];
		vectorRotations[69] = auxReferenceAxis[1];
		vectorRotations[70] = auxReferenceAxis[2];
		vectorRotations[71] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[72] = auxReferenceAxis[0];
		vectorRotations[73] = auxReferenceAxis[1];
		vectorRotations[74] = auxReferenceAxis[2];
		vectorRotations[75] = 4.0e0 * auxMath_._pi / 3.0e0;

		//c3 ---- 0 - 3 - 4
		auxReferenceAxis[0] = (mol0[0].x + mol0[3].x + mol0[4].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[3].y + mol0[4].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[3].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[76] = auxReferenceAxis[0];
		vectorRotations[77] = auxReferenceAxis[1];
		vectorRotations[78] = auxReferenceAxis[2];
		vectorRotations[79] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[80] = auxReferenceAxis[0];
		vectorRotations[81] = auxReferenceAxis[1];
		vectorRotations[82] = auxReferenceAxis[2];
		vectorRotations[83] = 4.0e0 * auxMath_._pi / 3.0e0;

		//c3 ---- 0 - 4 - 1
		auxReferenceAxis[0] = (mol0[0].x + mol0[4].x + mol0[1].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[4].y + mol0[1].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[4].z + mol0[1].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[84] = auxReferenceAxis[0];
		vectorRotations[85] = auxReferenceAxis[1];
		vectorRotations[86] = auxReferenceAxis[2];
		vectorRotations[87] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[88] = auxReferenceAxis[0];
		vectorRotations[89] = auxReferenceAxis[1];
		vectorRotations[90] = auxReferenceAxis[2];
		vectorRotations[91] = 4.0e0 * auxMath_._pi / 3.0e0;
		break;


	case 7:
		mol0.resize(7);
		mol0[0].x = 0.00000000;
		mol0[0].y = 0.00000000;
		mol0[0].z = 1.00000000;
		mol0[1].x = 0.97767167;
		mol0[1].y = 0.00000000;
		mol0[1].z = 0.21013831;
		mol0[2].x = 0.16977090;
		mol0[2].y = 0.96281864;
		mol0[2].z = 0.21013831;
		mol0[3].x = -0.91871085;
		mol0[3].y = 0.33438340;
		mol0[3].z = 0.21013831;
		mol0[4].x = -0.48883583;
		mol0[4].y = -0.84668850;
		mol0[4].z = 0.21013831;
		mol0[5].x = 0.36282725;
		mol0[5].y = -0.62843523;
		mol0[5].z = -0.68805926;
		mol0[6].x = -0.26010411;
		mol0[6].y = 0.45051354;
		mol0[6].z = -0.85403946;
		cutAngle = auxMath_._pi / 2.0e0;
		vectorRotations.resize(8);
		//c3 - 1
		vectorRotations[0] = mol0[5].x;
		vectorRotations[1] = mol0[5].y;
		vectorRotations[2] = mol0[5].z;
		vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
		//c3 - 2
		vectorRotations[4] = mol0[5].x;
		vectorRotations[5] = mol0[5].y;
		vectorRotations[6] = mol0[5].z;
		vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;
		break;



	case 8:
		mol0.resize(8);
		mol0[0].x = 0.000000;
		mol0[0].y = 0.000000;
		mol0[0].z = 1.00000000;
		mol0[1].x = 0.96528366;
		mol0[1].y = 0.000000;
		mol0[1].z = 0.26120388;
		mol0[2].x = -0.56545007;
		mol0[2].y = 0.78232905;
		mol0[2].z = 0.26120388;
		mol0[3].x = -0.88247541;
		mol0[3].y = -0.39116453;
		mol0[3].z = 0.26120387;
		mol0[4].x = 0.19991679;
		mol0[4].y = -0.94435471;
		mol0[4].z = 0.26120387;
		mol0[5].x = 0.39983358;
		mol0[5].y = 0.78232905;
		mol0[5].z = -0.47759225;
		mol0[6].x = -0.59975037;
		mol0[6].y = 0.16202565;
		mol0[6].z = -0.78361162;
		mol0[7].x = 0.48264183;
		mol0[7].y = -0.39116453;
		mol0[7].z = -0.78361162;
		cutAngle = auxMath_._pi / 2.0e0;
		vectorRotations.resize(28);
		auxReferenceAxis.resize(3);
		// 3 4 6 7
		auxReferenceAxis[0] = 0.25e0 *
			(mol0[3].x + mol0[4].x + mol0[6].x + mol0[7].x);
		auxReferenceAxis[1] = 0.25e0 *
			(mol0[3].y + mol0[4].y + mol0[6].y + mol0[7].y);
		auxReferenceAxis[2] = 0.25e0 *
			(mol0[3].z + mol0[4].z + mol0[6].z + mol0[7].z);
		auxMath_.normalize(auxReferenceAxis);
		//c4 - 1
		vectorRotations[0] = auxReferenceAxis[0];
		vectorRotations[1] = auxReferenceAxis[1];
		vectorRotations[2] = auxReferenceAxis[2];
		vectorRotations[3] = auxMath_._pi / 2.0e0;
		//c4 - 2 (c2)
		vectorRotations[4] = auxReferenceAxis[0];
		vectorRotations[5] = auxReferenceAxis[1];
		vectorRotations[6] = auxReferenceAxis[2];
		vectorRotations[7] = auxMath_._pi;
		//c4 - 3
		vectorRotations[8] = auxReferenceAxis[0];
		vectorRotations[9] = auxReferenceAxis[1];
		vectorRotations[10] = auxReferenceAxis[2];
		vectorRotations[11] = 3.0e0 * auxMath_._pi / 2.0e0;
		// c2 - 1 e 7
		auxReferenceAxis[0] = 0.5e0 *
			(mol0[1].x + mol0[7].x);
		auxReferenceAxis[1] = 0.5e0 *
			(mol0[1].y + mol0[7].y);
		auxReferenceAxis[2] = 0.5e0 *
			(mol0[1].z + mol0[7].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[12] = auxReferenceAxis[0];
		vectorRotations[13] = auxReferenceAxis[1];
		vectorRotations[14] = auxReferenceAxis[2];
		vectorRotations[15] = auxMath_._pi;
		// c2 - 7 e 5
		auxReferenceAxis[0] = 0.5e0 *
			(mol0[7].x + mol0[5].x);
		auxReferenceAxis[1] = 0.5e0 *
			(mol0[7].y + mol0[5].y);
		auxReferenceAxis[2] = 0.5e0 *
			(mol0[7].z + mol0[5].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[16] = auxReferenceAxis[0];
		vectorRotations[17] = auxReferenceAxis[1];
		vectorRotations[18] = auxReferenceAxis[2];
		vectorRotations[19] = auxMath_._pi;
		// c2 - 5 e 6
		auxReferenceAxis[0] = 0.5e0 *
			(mol0[5].x + mol0[6].x);
		auxReferenceAxis[1] = 0.5e0 *
			(mol0[5].y + mol0[6].y);
		auxReferenceAxis[2] = 0.5e0 *
			(mol0[5].z + mol0[6].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[20] = auxReferenceAxis[0];
		vectorRotations[21] = auxReferenceAxis[1];
		vectorRotations[22] = auxReferenceAxis[2];
		vectorRotations[23] = auxMath_._pi;
		// c2 - 6 e 2
		auxReferenceAxis[0] = 0.5e0 *
			(mol0[6].x + mol0[2].x);
		auxReferenceAxis[1] = 0.5e0 *
			(mol0[6].y + mol0[2].y);
		auxReferenceAxis[2] = 0.5e0 *
			(mol0[6].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[24] = auxReferenceAxis[0];
		vectorRotations[25] = auxReferenceAxis[1];
		vectorRotations[26] = auxReferenceAxis[2];
		vectorRotations[27] = auxMath_._pi;
		break;


	case 9:
		mol0.resize(9);
		mol0[0].x = 0.000000;
		mol0[0].y = 0.000000;
		mol0[0].z = 1.000000;
		mol0[1].x = -0.23570226;
		mol0[1].y = 0.91287093;
		mol0[1].z = 0.33333333;
		mol0[2].x = -0.94280904;
		mol0[2].y = 0.00000000;
		mol0[2].z = 0.33333333;
		mol0[3].x = 0.23570226;
		mol0[3].y = -0.91287093;
		mol0[3].z = 0.33333333;
		mol0[4].x = 0.94280904;
		mol0[4].y = 0.00000000;
		mol0[4].z = 0.33333333;
		mol0[5].x = 0.53033009;
		mol0[5].y = 0.68465320;
		mol0[5].z = -0.50000000;
		mol0[6].x = -0.53033009;
		mol0[6].y = -0.68465320;
		mol0[6].z = -0.50000000;
		mol0[7].x = -0.58925565;
		mol0[7].y = 0.45643546;
		mol0[7].z = -0.66666667;
		mol0[8].x = 0.58925565;
		mol0[8].y = -0.45643546;
		mol0[8].z = -0.66666667;
		cutAngle = auxMath_._pi / 2.0e0;
		vectorRotations.resize(20);
		//c2 - 1
		vectorRotations[0] = mol0[0].x;
		vectorRotations[1] = mol0[0].y;
		vectorRotations[2] = mol0[0].z;
		vectorRotations[3] = auxMath_._pi;
		//c2 - 2
		vectorRotations[4] = mol0[6].x;
		vectorRotations[5] = mol0[6].y;
		vectorRotations[6] = mol0[6].z;
		vectorRotations[7] = auxMath_._pi;
		//c2 - 3
		vectorRotations[8] = mol0[5].x;
		vectorRotations[9] = mol0[5].y;
		vectorRotations[10] = mol0[5].z;
		vectorRotations[11] = auxMath_._pi;
		// produto vetorial
		auxReferenceAxis = auxMath_.vectorProduct(
			mol0[0].x, mol0[0].y, mol0[0].z,
			mol0[6].x, mol0[6].y, mol0[6].z);
		auxMath_.normalize(auxReferenceAxis);
		//c3 - 1
		vectorRotations[12] = auxReferenceAxis[0];
		vectorRotations[13] = auxReferenceAxis[1];
		vectorRotations[14] = auxReferenceAxis[2];
		vectorRotations[15] = 2.0e0 * auxMath_._pi / 3.0e0;
		//c3 - 2
		vectorRotations[16] = auxReferenceAxis[0];
		vectorRotations[17] = auxReferenceAxis[1];
		vectorRotations[18] = auxReferenceAxis[2];
		vectorRotations[19] = 4.0e0 * auxMath_._pi / 3.0e0;
		break;


	case 10:
		mol0.resize(10);
		mol0[0].x = 0.000000;
		mol0[0].y = 0.000000;
		mol0[0].z = 1.000000;
		mol0[1].x = 0.91458473;
		mol0[1].y = 0.00000000;
		mol0[1].z = 0.40439433;
		mol0[2].x = 0.26335401;
		mol0[2].y = 0.87584810;
		mol0[2].z = 0.40439432;
		mol0[3].x = -0.76291954;
		mol0[3].y = 0.50439965;
		mol0[3].z = 0.40439433;
		mol0[4].x = -0.38787671;
		mol0[4].y = -0.82826136;
		mol0[4].z = 0.40439433;
		mol0[5].x = 0.52670802;
		mol0[5].y = -0.82826136;
		mol0[5].z = -0.19121135;
		mol0[6].x = -0.89351054;
		mol0[6].y = -0.25145533;
		mol0[6].z = -0.37203377;
		mol0[7].x = 0.67837321;
		mol0[7].y = 0.50439965;
		mol0[7].z = -0.53421979;
		mol0[8].x = -0.39464230;
		mol0[8].y = 0.69067213;
		mol0[8].z = -0.60599460;
		mol0[9].x = 0.02107419;
		mol0[9].y = -0.25145533;
		mol0[9].z = -0.96763944;
		cutAngle = auxMath_._pi / 2.0e0;
		vectorRotations.resize(4);
		//c2 - 1
		auxReferenceAxis.resize(3);
		auxReferenceAxis[0] = 0.5e0 * (mol0[7].x + mol0[3].x);
		auxReferenceAxis[1] = 0.5e0 * (mol0[7].y + mol0[3].y);
		auxReferenceAxis[2] = 0.5e0 * (mol0[7].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[0] = auxReferenceAxis[0];
		vectorRotations[1] = auxReferenceAxis[1];
		vectorRotations[2] = auxReferenceAxis[2];
		vectorRotations[3] = auxMath_._pi;
		break;


	case 11:
		mol0.resize(11);
		mol0[0].x = 0.000000;
		mol0[0].y = 0.000000;
		mol0[0].z = 1.000000;
		mol0[1].x = 0.89442719;
		mol0[1].y = 0.00000000;
		mol0[1].z = 0.44721360;
		mol0[2].x = 0.27639320;
		mol0[2].y = 0.85065081;
		mol0[2].z = 0.44721360;
		mol0[3].x = -0.72360680;
		mol0[3].y = 0.52573111;
		mol0[3].z = 0.44721360;
		mol0[4].x = -0.72360680;
		mol0[4].y = -0.52573111;
		mol0[4].z = 0.44721360;
		mol0[5].x = 0.27639320;
		mol0[5].y = -0.85065081;
		mol0[5].z = 0.44721360;
		mol0[6].x = 0.72360680;
		mol0[6].y = 0.52573111;
		mol0[6].z = -0.44721360;
		mol0[7].x = -0.27639320;
		mol0[7].y = 0.85065081;
		mol0[7].z = -0.44721360;
		mol0[8].x = -0.89442719;
		mol0[8].y = 0.00000000;
		mol0[8].z = -0.44721360;
		mol0[9].x = -0.27639320;
		mol0[9].y = -0.85065081;
		mol0[9].z = -0.44721360;
		mol0[10].x = 0.00000000;
		mol0[10].y = 0.00000000;
		mol0[10].z = -1.00000000;
		cutAngle = auxMath_._pi / 2.0e0;
		vectorRotations.resize(16);
		//add rotations
		//c5 - 1
		vectorRotations[0] = mol0[3].x;
		vectorRotations[1] = mol0[3].y;
		vectorRotations[2] = mol0[3].z;
		vectorRotations[3] = (2.0e0 * auxMath_._pi / 5.0e0);
		//c5 - 2
		vectorRotations[4] = mol0[3].x;
		vectorRotations[5] = mol0[3].y;
		vectorRotations[6] = mol0[3].z;
		vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
		//c5 - 3
		vectorRotations[8] = mol0[3].x;
		vectorRotations[9] = mol0[3].y;
		vectorRotations[10] = mol0[3].z;
		vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
		//c5 - 4
		vectorRotations[12] = mol0[3].x;
		vectorRotations[13] = mol0[3].y;
		vectorRotations[14] = mol0[3].z;
		vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
		break;


	case 12:
		mol0.resize(12);
		mol0[0].x = 0.00000000;
		mol0[0].y = 0.00000000;
		mol0[0].z = 1.00000000;
		mol0[1].x = 0.89442719;
		mol0[1].y = 0.00000000;
		mol0[1].z = 0.44721360;
		mol0[2].x = 0.27639320;
		mol0[2].y = 0.85065081;
		mol0[2].z = 0.44721360;
		mol0[3].x = -0.72360680;
		mol0[3].y = 0.52573111;
		mol0[3].z = 0.44721360;
		mol0[4].x = -0.72360680;
		mol0[4].y = -0.52573111;
		mol0[4].z = 0.44721360;
		mol0[5].x = 0.27639320;
		mol0[5].y = -0.85065081;
		mol0[5].z = 0.44721360;
		mol0[6].x = 0.72360680;
		mol0[6].y = 0.52573111;
		mol0[6].z = -0.44721360;
		mol0[7].x = -0.27639320;
		mol0[7].y = 0.85065081;
		mol0[7].z = -0.44721360;
		mol0[8].x = -0.89442719;
		mol0[8].y = 0.00000000;
		mol0[8].z = -0.44721360;
		mol0[9].x = -0.27639320;
		mol0[9].y = -0.85065081;
		mol0[9].z = -0.44721360;
		mol0[10].x = 0.72360680;
		mol0[10].y = -0.52573111;
		mol0[10].z = -0.44721360;
		mol0[11].x = 0.00000000;
		mol0[11].y = 0.00000000;
		mol0[11].z = -1.00000000;
		cutAngle = auxMath_._pi / 2.0e0;
		vectorRotations.resize(236);
		// 4 c5 para cada um dos 6 pontos principais (12 c5 e 12 c52)
		// 2 c3 pra cada face unica (20 c3)
		// c2 para cada aresta (10 no topo e 5 ao redor = 15)
		//add rotations
		// 0 por cima
		// c5 - 0
		vectorRotations[0] = mol0[0].x;
		vectorRotations[1] = mol0[0].y;
		vectorRotations[2] = mol0[0].z;
		vectorRotations[3] = (2.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[4] = mol0[0].x;
		vectorRotations[5] = mol0[0].y;
		vectorRotations[6] = mol0[0].z;
		vectorRotations[7] = (4.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[8] = mol0[0].x;
		vectorRotations[9] = mol0[0].y;
		vectorRotations[10] = mol0[0].z;
		vectorRotations[11] = (6.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[12] = mol0[0].x;
		vectorRotations[13] = mol0[0].y;
		vectorRotations[14] = mol0[0].z;
		vectorRotations[15] = (8.0e0 * auxMath_._pi / 5.0e0);

		// c5 - 1
		vectorRotations[16] = mol0[1].x;
		vectorRotations[17] = mol0[1].y;
		vectorRotations[18] = mol0[1].z;
		vectorRotations[19] = (2.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[20] = mol0[1].x;
		vectorRotations[21] = mol0[1].y;
		vectorRotations[22] = mol0[1].z;
		vectorRotations[23] = (4.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[24] = mol0[1].x;
		vectorRotations[25] = mol0[1].y;
		vectorRotations[26] = mol0[1].z;
		vectorRotations[27] = (6.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[28] = mol0[1].x;
		vectorRotations[29] = mol0[1].y;
		vectorRotations[30] = mol0[1].z;
		vectorRotations[31] = (8.0e0 * auxMath_._pi / 5.0e0);

		// c5 - 2
		vectorRotations[32] = mol0[2].x;
		vectorRotations[33] = mol0[2].y;
		vectorRotations[34] = mol0[2].z;
		vectorRotations[35] = (2.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[36] = mol0[2].x;
		vectorRotations[37] = mol0[2].y;
		vectorRotations[38] = mol0[2].z;
		vectorRotations[39] = (4.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[40] = mol0[2].x;
		vectorRotations[41] = mol0[2].y;
		vectorRotations[42] = mol0[2].z;
		vectorRotations[43] = (6.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[44] = mol0[2].x;
		vectorRotations[45] = mol0[2].y;
		vectorRotations[46] = mol0[2].z;
		vectorRotations[47] = (8.0e0 * auxMath_._pi / 5.0e0);

		// c5 - 3
		vectorRotations[48] = mol0[3].x;
		vectorRotations[49] = mol0[3].y;
		vectorRotations[50] = mol0[3].z;
		vectorRotations[51] = (2.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[52] = mol0[3].x;
		vectorRotations[53] = mol0[3].y;
		vectorRotations[54] = mol0[3].z;
		vectorRotations[55] = (4.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[56] = mol0[3].x;
		vectorRotations[57] = mol0[3].y;
		vectorRotations[58] = mol0[3].z;
		vectorRotations[59] = (6.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[60] = mol0[3].x;
		vectorRotations[61] = mol0[3].y;
		vectorRotations[62] = mol0[3].z;
		vectorRotations[63] = (8.0e0 * auxMath_._pi / 5.0e0);

		// c5 - 4
		vectorRotations[64] = mol0[4].x;
		vectorRotations[65] = mol0[4].y;
		vectorRotations[66] = mol0[4].z;
		vectorRotations[67] = (2.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[68] = mol0[4].x;
		vectorRotations[69] = mol0[4].y;
		vectorRotations[70] = mol0[4].z;
		vectorRotations[71] = (4.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[72] = mol0[4].x;
		vectorRotations[73] = mol0[4].y;
		vectorRotations[74] = mol0[4].z;
		vectorRotations[75] = (6.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[76] = mol0[4].x;
		vectorRotations[77] = mol0[4].y;
		vectorRotations[78] = mol0[4].z;
		vectorRotations[79] = (8.0e0 * auxMath_._pi / 5.0e0);

		// c5 - 5
		vectorRotations[80] = mol0[5].x;
		vectorRotations[81] = mol0[5].y;
		vectorRotations[82] = mol0[5].z;
		vectorRotations[83] = (2.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[84] = mol0[5].x;
		vectorRotations[85] = mol0[5].y;
		vectorRotations[86] = mol0[5].z;
		vectorRotations[87] = (4.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[88] = mol0[5].x;
		vectorRotations[89] = mol0[5].y;
		vectorRotations[90] = mol0[5].z;
		vectorRotations[91] = (6.0e0 * auxMath_._pi / 5.0e0);
		vectorRotations[92] = mol0[5].x;
		vectorRotations[93] = mol0[5].y;
		vectorRotations[94] = mol0[5].z;
		vectorRotations[95] = (8.0e0 * auxMath_._pi / 5.0e0);

		auxReferenceAxis.resize(3);
		// c3 pra cada face
		// 0-1-2
		auxReferenceAxis[0] = (mol0[0].x + mol0[1].x + mol0[2].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[1].y + mol0[2].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[1].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[96] = auxReferenceAxis[0];
		vectorRotations[97] = auxReferenceAxis[1];
		vectorRotations[98] = auxReferenceAxis[2];
		vectorRotations[99] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[100] = auxReferenceAxis[0];
		vectorRotations[101] = auxReferenceAxis[1];
		vectorRotations[102] = auxReferenceAxis[2];
		vectorRotations[103] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 0-2-3
		auxReferenceAxis[0] = (mol0[0].x + mol0[2].x + mol0[3].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[2].y + mol0[3].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[2].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[104] = auxReferenceAxis[0];
		vectorRotations[105] = auxReferenceAxis[1];
		vectorRotations[106] = auxReferenceAxis[2];
		vectorRotations[107] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[108] = auxReferenceAxis[0];
		vectorRotations[109] = auxReferenceAxis[1];
		vectorRotations[110] = auxReferenceAxis[2];
		vectorRotations[111] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 0-3-4
		auxReferenceAxis[0] = (mol0[0].x + mol0[3].x + mol0[4].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[3].y + mol0[4].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[3].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[112] = auxReferenceAxis[0];
		vectorRotations[113] = auxReferenceAxis[1];
		vectorRotations[114] = auxReferenceAxis[2];
		vectorRotations[115] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[116] = auxReferenceAxis[0];
		vectorRotations[117] = auxReferenceAxis[1];
		vectorRotations[118] = auxReferenceAxis[2];
		vectorRotations[119] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 0-4-5
		auxReferenceAxis[0] = (mol0[0].x + mol0[4].x + mol0[5].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[4].y + mol0[5].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[4].z + mol0[5].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[120] = auxReferenceAxis[0];
		vectorRotations[121] = auxReferenceAxis[1];
		vectorRotations[122] = auxReferenceAxis[2];
		vectorRotations[123] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[124] = auxReferenceAxis[0];
		vectorRotations[125] = auxReferenceAxis[1];
		vectorRotations[126] = auxReferenceAxis[2];
		vectorRotations[127] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 0-5-1
		auxReferenceAxis[0] = (mol0[0].x + mol0[5].x + mol0[1].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[5].y + mol0[1].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[5].z + mol0[1].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[128] = auxReferenceAxis[0];
		vectorRotations[129] = auxReferenceAxis[1];
		vectorRotations[130] = auxReferenceAxis[2];
		vectorRotations[131] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[132] = auxReferenceAxis[0];
		vectorRotations[133] = auxReferenceAxis[1];
		vectorRotations[134] = auxReferenceAxis[2];
		vectorRotations[135] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 3-8-4
		auxReferenceAxis[0] = (mol0[3].x + mol0[8].x + mol0[4].x);
		auxReferenceAxis[1] = (mol0[3].y + mol0[8].y + mol0[4].y);
		auxReferenceAxis[2] = (mol0[3].z + mol0[8].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[136] = auxReferenceAxis[0];
		vectorRotations[137] = auxReferenceAxis[1];
		vectorRotations[138] = auxReferenceAxis[2];
		vectorRotations[139] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[140] = auxReferenceAxis[0];
		vectorRotations[141] = auxReferenceAxis[1];
		vectorRotations[142] = auxReferenceAxis[2];
		vectorRotations[143] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 8-4-9
		auxReferenceAxis[0] = (mol0[8].x + mol0[4].x + mol0[9].x);
		auxReferenceAxis[1] = (mol0[8].y + mol0[4].y + mol0[9].y);
		auxReferenceAxis[2] = (mol0[8].z + mol0[4].z + mol0[9].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[144] = auxReferenceAxis[0];
		vectorRotations[145] = auxReferenceAxis[1];
		vectorRotations[146] = auxReferenceAxis[2];
		vectorRotations[147] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[148] = auxReferenceAxis[0];
		vectorRotations[149] = auxReferenceAxis[1];
		vectorRotations[150] = auxReferenceAxis[2];
		vectorRotations[151] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 4-9-5
		auxReferenceAxis[0] = (mol0[4].x + mol0[9].x + mol0[5].x);
		auxReferenceAxis[1] = (mol0[4].y + mol0[9].y + mol0[5].y);
		auxReferenceAxis[2] = (mol0[4].z + mol0[9].z + mol0[5].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[152] = auxReferenceAxis[0];
		vectorRotations[153] = auxReferenceAxis[1];
		vectorRotations[154] = auxReferenceAxis[2];
		vectorRotations[155] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[156] = auxReferenceAxis[0];
		vectorRotations[157] = auxReferenceAxis[1];
		vectorRotations[158] = auxReferenceAxis[2];
		vectorRotations[159] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 9-5-10
		auxReferenceAxis[0] = (mol0[9].x + mol0[5].x + mol0[10].x);
		auxReferenceAxis[1] = (mol0[9].y + mol0[5].y + mol0[10].y);
		auxReferenceAxis[2] = (mol0[9].z + mol0[5].z + mol0[10].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[160] = auxReferenceAxis[0];
		vectorRotations[161] = auxReferenceAxis[1];
		vectorRotations[162] = auxReferenceAxis[2];
		vectorRotations[163] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[164] = auxReferenceAxis[0];
		vectorRotations[165] = auxReferenceAxis[1];
		vectorRotations[166] = auxReferenceAxis[2];
		vectorRotations[167] = 4.0e0 * auxMath_._pi / 3.0e0;

		// 5-10-1
		auxReferenceAxis[0] = (mol0[10].x + mol0[5].x + mol0[1].x);
		auxReferenceAxis[1] = (mol0[10].y + mol0[5].y + mol0[1].y);
		auxReferenceAxis[2] = (mol0[10].z + mol0[5].z + mol0[1].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[168] = auxReferenceAxis[0];
		vectorRotations[169] = auxReferenceAxis[1];
		vectorRotations[170] = auxReferenceAxis[2];
		vectorRotations[171] = 2.0e0 * auxMath_._pi / 3.0e0;
		vectorRotations[172] = auxReferenceAxis[0];
		vectorRotations[173] = auxReferenceAxis[1];
		vectorRotations[174] = auxReferenceAxis[2];
		vectorRotations[175] = 4.0e0 * auxMath_._pi / 3.0e0;

		// c2 nas arestas
		// 0-1
		auxReferenceAxis[0] = (mol0[0].x + mol0[1].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[1].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[1].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[176] = auxReferenceAxis[0];
		vectorRotations[177] = auxReferenceAxis[1];
		vectorRotations[178] = auxReferenceAxis[2];
		vectorRotations[179] = auxMath_._pi;

		// 0-2
		auxReferenceAxis[0] = (mol0[0].x + mol0[2].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[2].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[180] = auxReferenceAxis[0];
		vectorRotations[181] = auxReferenceAxis[1];
		vectorRotations[182] = auxReferenceAxis[2];
		vectorRotations[183] = auxMath_._pi;

		// 0-3
		auxReferenceAxis[0] = (mol0[0].x + mol0[3].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[3].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[184] = auxReferenceAxis[0];
		vectorRotations[185] = auxReferenceAxis[1];
		vectorRotations[186] = auxReferenceAxis[2];
		vectorRotations[187] = auxMath_._pi;

		// 0-4
		auxReferenceAxis[0] = (mol0[0].x + mol0[4].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[4].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[188] = auxReferenceAxis[0];
		vectorRotations[189] = auxReferenceAxis[1];
		vectorRotations[190] = auxReferenceAxis[2];
		vectorRotations[191] = auxMath_._pi;

		// 0-5
		auxReferenceAxis[0] = (mol0[0].x + mol0[5].x);
		auxReferenceAxis[1] = (mol0[0].y + mol0[5].y);
		auxReferenceAxis[2] = (mol0[0].z + mol0[5].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[192] = auxReferenceAxis[0];
		vectorRotations[193] = auxReferenceAxis[1];
		vectorRotations[194] = auxReferenceAxis[2];
		vectorRotations[195] = auxMath_._pi;

		// 1-2
		auxReferenceAxis[0] = (mol0[1].x + mol0[2].x);
		auxReferenceAxis[1] = (mol0[1].y + mol0[2].y);
		auxReferenceAxis[2] = (mol0[1].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[196] = auxReferenceAxis[0];
		vectorRotations[197] = auxReferenceAxis[1];
		vectorRotations[198] = auxReferenceAxis[2];
		vectorRotations[199] = auxMath_._pi;

		// 2-3
		auxReferenceAxis[0] = (mol0[2].x + mol0[3].x);
		auxReferenceAxis[1] = (mol0[2].y + mol0[3].y);
		auxReferenceAxis[2] = (mol0[2].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[200] = auxReferenceAxis[0];
		vectorRotations[201] = auxReferenceAxis[1];
		vectorRotations[202] = auxReferenceAxis[2];
		vectorRotations[203] = auxMath_._pi;

		// 3-4
		auxReferenceAxis[0] = (mol0[3].x + mol0[4].x);
		auxReferenceAxis[1] = (mol0[3].y + mol0[4].y);
		auxReferenceAxis[2] = (mol0[3].z + mol0[4].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[204] = auxReferenceAxis[0];
		vectorRotations[205] = auxReferenceAxis[1];
		vectorRotations[206] = auxReferenceAxis[2];
		vectorRotations[207] = auxMath_._pi;

		// 4-5
		auxReferenceAxis[0] = (mol0[4].x + mol0[5].x);
		auxReferenceAxis[1] = (mol0[4].y + mol0[5].y);
		auxReferenceAxis[2] = (mol0[4].z + mol0[5].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[208] = auxReferenceAxis[0];
		vectorRotations[209] = auxReferenceAxis[1];
		vectorRotations[210] = auxReferenceAxis[2];
		vectorRotations[211] = auxMath_._pi;

		// 5-1
		auxReferenceAxis[0] = (mol0[5].x + mol0[1].x);
		auxReferenceAxis[1] = (mol0[5].y + mol0[1].y);
		auxReferenceAxis[2] = (mol0[5].z + mol0[1].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[212] = auxReferenceAxis[0];
		vectorRotations[213] = auxReferenceAxis[1];
		vectorRotations[214] = auxReferenceAxis[2];
		vectorRotations[215] = auxMath_._pi;

		// 1-6
		auxReferenceAxis[0] = (mol0[1].x + mol0[6].x);
		auxReferenceAxis[1] = (mol0[1].y + mol0[6].y);
		auxReferenceAxis[2] = (mol0[1].z + mol0[6].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[216] = auxReferenceAxis[0];
		vectorRotations[217] = auxReferenceAxis[1];
		vectorRotations[218] = auxReferenceAxis[2];
		vectorRotations[219] = auxMath_._pi;

		// 6-2
		auxReferenceAxis[0] = (mol0[6].x + mol0[2].x);
		auxReferenceAxis[1] = (mol0[6].y + mol0[2].y);
		auxReferenceAxis[2] = (mol0[6].z + mol0[2].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[220] = auxReferenceAxis[0];
		vectorRotations[221] = auxReferenceAxis[1];
		vectorRotations[222] = auxReferenceAxis[2];
		vectorRotations[223] = auxMath_._pi;

		// 2-7
		auxReferenceAxis[0] = (mol0[2].x + mol0[7].x);
		auxReferenceAxis[1] = (mol0[2].y + mol0[7].y);
		auxReferenceAxis[2] = (mol0[2].z + mol0[7].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[224] = auxReferenceAxis[0];
		vectorRotations[225] = auxReferenceAxis[1];
		vectorRotations[226] = auxReferenceAxis[2];
		vectorRotations[227] = auxMath_._pi;

		// 7-3
		auxReferenceAxis[0] = (mol0[7].x + mol0[3].x);
		auxReferenceAxis[1] = (mol0[7].y + mol0[3].y);
		auxReferenceAxis[2] = (mol0[7].z + mol0[3].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[228] = auxReferenceAxis[0];
		vectorRotations[229] = auxReferenceAxis[1];
		vectorRotations[230] = auxReferenceAxis[2];
		vectorRotations[231] = auxMath_._pi;

		// 3-8
		auxReferenceAxis[0] = (mol0[3].x + mol0[8].x);
		auxReferenceAxis[1] = (mol0[3].y + mol0[8].y);
		auxReferenceAxis[2] = (mol0[3].z + mol0[8].z);
		auxMath_.normalize(auxReferenceAxis);
		vectorRotations[232] = auxReferenceAxis[0];
		vectorRotations[233] = auxReferenceAxis[1];
		vectorRotations[234] = auxReferenceAxis[2];
		vectorRotations[235] = auxMath_._pi;
		break;

	default:
		cout << "CauchyIndex::setSystem - system not found" << endl;
		exit(1);
		break;
	}
	setAllRotations(vectorRotations);
}















/* TESTANDO IO
for (size_t i = 0; i < allPermutations.size(); i++)
{
writeCauchyRotations("allPermut.txt", allPermutations[i].rotPermutations);
}

ifstream ifOpen_("allPermut.txt");
for (size_t i = 0; i < 120; i++)
{
vector<int> permutTeste = readCauchyNotations(ifOpen_);
printCauchyNotation(permutTeste);
}
ifOpen_.close();
*/

/*
bidentateMap[0][1] = 1;
bidentateMap[0][2] = 7;
bidentateMap[0][3] = 4;
bidentateMap[0][4] = 5;
bidentateMap[1][2] = -1;
bidentateMap[1][3] = 2;
bidentateMap[1][4] = 3;
bidentateMap[2][3] = 8;
bidentateMap[2][4] = 9;
bidentateMap[3][4] = 6;
*/



/* CODIGO Q DEVE SERVIR PRA ALGUMA COISA
vector<int> permutation(5);
vector<string> atoms(5);
for (size_t i = 0; i < 5; i++)
{
permutation[i] = i;
atoms[i] = "H ";
}
atoms[0] = "C ";
atoms[1] = "C ";
vector<int> bidentateAtomsChosen(2);
bidentateAtomsChosen[0] = 0;
bidentateAtomsChosen[1] = 1;

ci_.rotationTest(atoms, bidentateAtomsChosen);

ofstream of_("printCauchy.xyz");
ci_.printMolecule(permutation, atoms, bidentateAtomsChosen,of_);
of_.close();
*/
