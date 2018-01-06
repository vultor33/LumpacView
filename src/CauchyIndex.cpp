#include "CauchyIndex.h"

//#define UNIX

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
#include "IdentifyIsomers.h"
#include "Geometries.h"
#include "IsomersToMol.h"

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
	atomLabels[0] = "Ca";
	atomLabels[1] = "O";
	atomLabels[2] = "H";
	atomLabels[3] = "Mg";
	atomLabels[4] = "Be";
	atomLabels[5] = "I";
	atomLabels[6] = "Cl";
	atomLabels[7] = "Br";
	atomLabels[8] = "Na";
	atomLabels[9] = "He";
	atomLabels[10] = "N";
	atomLabels[11] = "C";


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


// nao tem o que fazer - toda vez que eu deletar manda uma  linha
// em um arquivo -> k
void CauchyIndex::generateAllIndependentIsomersWithFlagEnantiomers(
	string blockFileName,
	string flagsFile,
	string code)
{

	vector<int> atomTypes;
	vector<int> bidentateAtomChosen;
	readAtomTypesAndBidentateChosenFile(flagsFile, atomTypes, bidentateAtomChosen);

	systemSize = mol0.size();
	int size = factorial(systemSize);
	nRotations = allRotationTransforms.size() + 1;
	vector< vector<int> > allPermutations;
	vector< vector<int> > blockPermutations;
	vector<int> weights;
	ifstream openedFile_(blockFileName.c_str());
	int compareResult;
	bool equal;
	int k = 0;
	int kDupe = 0;
	ofstream of_((code + "-" + blockFileName).c_str());
	while (!openedFile_.eof())
	{
		vector<int> auxBlock = readCauchyNotations(openedFile_);
		if (auxBlock.size() == 0)
		{
			if (kDupe != 0)
			{
				of_ << endl;
				kDupe = 0;
			}
			continue;
		}
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
			{
				if(compareResult == 0)
					weights[i]++;
				break;
			}
		}
		if (!equal)
		{
			allPermutations.push_back(auxBlock);
			weights.push_back(0);
			for (size_t j = 0; j < systemSize; j++)
			{
				of_ << allPermutations[k][j] << "  ";
			}
			of_ << endl;
			k++;
			kDupe++;
		}
	}
	of_.close();

	string blockNum = blockFileName.substr(20,blockFileName.size());
	stringstream convBlockNum;
	convBlockNum << blockNum;
	int iBlock;
	convBlockNum >> iBlock;
	stringstream numberBackToString;
	numberBackToString << iBlock;
	if(iBlock < 10)
		blockNum = "000";
	else if(iBlock < 100)
		blockNum = "00";
	else if(iBlock < 1000)
		blockNum = "0";
	blockNum = blockNum + numberBackToString.str();

	ofstream ofCsv_((code + "-" + blockNum + ".csv").c_str());
	ifstream fileEnant_((code + "-" + blockFileName).c_str());
	int i = 0;
	while(true)
	{
		vector<int> auxEnant = readCauchyNotations(fileEnant_);
		if (fileEnant_.eof())
		{
			if(i != 0)
				ofCsv_ << endl;
			fileEnant_.close();
			break;
		}
		if (auxEnant.size() == 0)
		{
			ofCsv_ << endl;
			continue;
		}
		ofCsv_ << weights[i] << " ; ";
		for (size_t j = 0; j < systemSize; j++)
		{
			ofCsv_ << auxEnant[j] << "  ";
		}
		ofCsv_ << endl;
		i++;
		if (i == weights.size())
		{
			if(i != 0)
				ofCsv_ << endl;
			break;
		}
	}
	ofCsv_.close();
	fileEnant_.close();
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

	if(permutationIsomer1 == permutationIsomer2)
		return true;

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

	int totalSize = factorial(mol0.size());
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
                mol0.size(),
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
		if(haveToBreak)
			break;
	}


	ofstream isomers1_("independent-isomers-1");
        ifstream openedFile_(rawIsomersFile.c_str());
	for(int i = 1; i <= bigBlockSize; i++)
	{
                vector<int> permutation = readCauchyNotations(openedFile_);
		if(permutation.size() == 0)
		{
			i--;
			isomers1_ << endl;
			continue;
		}
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


void CauchyIndex::doBlockDeletionFlagsEnantiomers(
	string skeletonFile,
	string flagsFile,
	int kInit,
	int kFinal,
	int ramInit,
	int ramFinal)
{
	int systemSize = mol0.size();
	ifstream openedFile_(skeletonFile.c_str());
	vector< vector<int> > ramBlock = readCauchyNotationsRAMBlockEnantiomers(openedFile_, ramInit, ramFinal);

	openedFile_.close();
	vector<int> atomTypes;
	vector<int> bidentateChosen;
	readAtomTypesAndBidentateChosenFile(flagsFile, atomTypes, bidentateChosen);
	openedFile_.open(skeletonFile.c_str());
	ofstream weightFile_;
	stringstream weightIndex;
	weightIndex << ramInit << "to" << ramFinal;
	string weightFileName = flagsFile + "---weights---" + weightIndex.str() + ".txt";
	int k = 1;
	do
	{
		vector<int> permutation = readCauchyNotations(openedFile_);
		if (permutation.size() == 0)
			continue;
		if (k >= kInit && k <= kFinal) // from block use this to delete
		{
			int compare;
			int kRam = 0;
			while (kRam < ramBlock.size())
			{
				if(ramBlock[kRam][0] == -1)
				{
					kRam++;
					continue;
				}
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
					if(compare == 0)
					{
						weightFile_.open((weightFileName).c_str(),  std::ofstream::out | std::ofstream::app);
						weightFile_ << k << endl;
						weightFile_.close();
					}
				}
				kRam++;
			}
			for (size_t j = 0; j < allRotationTransforms.size(); j++)
			{
				vector<int> auxPerm = applyRotation(permutation, j);
				kRam = 0;
				while(kRam < ramBlock.size())
				{
					if(ramBlock[kRam][0] == -1)
					{
						kRam++;
						continue;
					}
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
						if(compare == 0)
						{
							weightFile_.open((weightFileName).c_str(),  std::ofstream::out | std::ofstream::app);
							weightFile_ << k << endl;
							weightFile_.close();
						}
					}
					kRam++;
				}
			}
		}
		k++;
		if(k > kFinal)
			break;
	} while (true);
	openedFile_.close();

	stringstream convert;
	convert << "block-" << ramInit << "-to-" << ramFinal << ".txt";
	string blockDelName = convert.str();
	ofstream of_(blockDelName.c_str());
	for (size_t i = 0; i < ramBlock.size(); i++)
	{
		if(ramBlock[i][0] == -1)
		{
			of_ << endl;
			continue;
		}
		for (size_t j = 0; j < systemSize; j++)
		{
			of_ << ramBlock[i][j] << "  ";
		}
		of_ << endl;
	}
	of_.close();
}




void CauchyIndex::createEnantiomersFiles(
	int nSystem,
	int nProc,
	int iMax,
	std::string skeletonFile,
	std::string workingDirectory,
	std::string machineType)
{
	for (int j = 1; j <= nProc; j++)
	{
		int i = j;
		stringstream convert;
		convert << j;
		string name_;
		if(machineType == "pc")
			name_ = convert.str() + ".x";
		else
			name_ = convert.str() + ".srm";
		ofstream roda_(name_.c_str());
		if (machineType == "pc")
			roda_ << "#!/bin/bash" << endl;
		else
		{
			roda_ << "#!/bin/bash" << endl;
			roda_ << "#SBATCH -n 1" << endl;
			roda_ << "#SBATCH --hint=nomultithread" << endl;
			roda_ << "RODADIR=" << workingDirectory << endl;
			roda_ << "cd $RODADIR" << endl;
		}
		while (i < iMax)
		{
			roda_ << "./lumpacview.exe enantiomerDeletion "
				<< nSystem << "  "
				<< i << "  "
				<< iMax << "  "
				<< convert.str() << "  "
				<< skeletonFile << endl;
			i += nProc;
		}
		roda_.close();
		system(("chmod u+x " + name_).c_str());
	}
}

void CauchyIndex::temporario()
{
	//"0 2 3 4 5 9 1 7 8 6"
	vector<int> isomerReflection(10);
	isomerReflection[0] = 0;
	isomerReflection[1] = 2;
	isomerReflection[2] = 3;
	isomerReflection[3] = 4;
	isomerReflection[4] = 5;
	isomerReflection[5] = 9;
	isomerReflection[6] = 1;
	isomerReflection[7] = 7;
	isomerReflection[8] = 8;
	isomerReflection[9] = 6;

	vector<int> distortedI = applyZAxisReflection(isomerReflection);

	printCauchyNotation(distortedI);

	cout << "ROTACOES" << endl;
	size_t nRotations = allRotationTransforms.size();
	for (size_t j = 0; j < nRotations; j++)
	{
		vector<int> auxPerm = applyRotation(isomerReflection, j);
		printCauchyNotation(auxPerm);
	}




}


void CauchyIndex::enantiomersOrdering()
{
	ifstream openedFile_("skeleton-isomers.txt");

	vector< vector<int> > allIsomers;
	while (true)
	{
		vector<int> isomer = readCauchyNotations(openedFile_);
		if (isomer.size() == 0)
			break;

		allIsomers.push_back(isomer);
	}
	openedFile_.close();

	vector<bool> taken(allIsomers.size());
	for (size_t i = 0; i < allIsomers.size(); i++)
		taken[i] = false;

	string fileName = "enatiomers-selection.txt";
	ofstream enantiomersOrder_(fileName.c_str());
	for (size_t i = 0; i < allIsomers.size(); i++)
	{
		if (taken[i])
			continue;
		// reflection
		vector<int> distortedI = applyZAxisReflection(allIsomers[i]);
		for (size_t j = i; j < allIsomers.size(); j++)
		{
			if (taken[j])
				continue;

			bool equal = compareTwoIsomers(distortedI, allIsomers[j]);

			if ((i == j) && equal)
			{
				for (size_t k = 0; k < allIsomers[i].size(); k++)
				{
					enantiomersOrder_ << allIsomers[i][k] << " ";
				}
				enantiomersOrder_ << endl << endl;
				taken[i] = true;
				break;
			}
			else if (equal)
			{
				for (size_t k = 0; k < allIsomers[i].size(); k++)
				{
					enantiomersOrder_ << allIsomers[i][k] << " ";
				}
				enantiomersOrder_ << endl;
				for (size_t k = 0; k < allIsomers[j].size(); k++)
				{
					enantiomersOrder_ << allIsomers[j][k] << " ";
				}
				enantiomersOrder_ << endl << endl;
				taken[i] = true;
				taken[j] = true;
				break;
			}
		}
	}
}


void CauchyIndex::enantiomersOrderingBlock(
	int iPermut,
	int kFinal,
	int nProc,
	string fileName)
{
	ifstream openedFile_(fileName.c_str());
	vector<int> distortedI;
	vector<int> isomer;
	vector<int> permutation;
	stringstream convert;
	convert << nProc;
	string numberProc = convert.str();
	int kInit = iPermut + 1;
	int k = 1;
	while (true)
	{
		isomer = readCauchyNotations(openedFile_);
		if (k == iPermut)
		{
			permutation = isomer;
			distortedI = applyZAxisReflection(permutation);
		}
		if (k >= kInit && k <= kFinal)
		{
			bool equal = compareTwoIsomers(distortedI, isomer);
			if (equal)
			{
				ofstream enantiomersOrder_;
				enantiomersOrder_.open(("enantiomers---" + numberProc + "-" + fileName).c_str(), std::ofstream::out | std::ofstream::app);
				for (size_t k = 0; k < permutation.size(); k++)
				{
					enantiomersOrder_ << permutation[k] << " ";
				}
				enantiomersOrder_ << endl;
				for (size_t k = 0; k < isomer.size(); k++)
				{
					enantiomersOrder_ << isomer[k] << " ";
				}
				enantiomersOrder_ << endl << endl;
				enantiomersOrder_.close();
				break;
			}
		}
		if (k > kFinal)
			break;
		k++;
	}
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
	ofstream of_("printTestRotations.xyz"); // reflection
	for (size_t i = 0; i < allRotationTransforms.size(); i++)
	{
		printMolecule(permutation, atoms, bidentateAtomsChosen, of_);
		vector<int> auxPerm = applyRotation(permutation, i);
		printMolecule(auxPerm, atoms, bidentateAtomsChosen, of_);
		printCauchyNotation(auxPerm);
		cout << endl;
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


vector<int> CauchyIndex::applyZAxisReflection(vector<int> permutation)
{
	vector<int> permutReflected = permutation;

	for (size_t i = 0; i < permutation.size(); i++)
	{
		permutReflected[reflectionOperation[i]] = permutation[i];
	}
	return permutReflected;
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
	printFile_ << setprecision(12);
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

void CauchyIndex::printMoleculeMolFormat(
	vector<int> & permutation,
	vector<int> & atomTypes,
	const vector<int> & bidentateAtomsChosen,
	ofstream & printFile_)
{
	double scale = 3.0e0;
	vector<string> atoms(permutation.size());
	for (size_t i = 0; i < atomTypes.size(); i++)
		atoms[i] = atomLabels[atomTypes[i]];

	if (permutation.size() == 0)
	{
		permutation.resize(atoms.size());
		for (size_t i = 0; i < atoms.size(); i++)
			permutation[i] = i;
	}
	//defining bidentate bonds
	vector<int> bidentateAtomsChosenRotated = applyPermutationCoordinates(permutation, atoms, bidentateAtomsChosen);
	int nAtoms = mol0.size();
	int nBonds = nAtoms + (bidentateAtomsChosenRotated.size() / 2);
	printFile_ <<"#   generated by: ComplexBuild" << endl << endl << endl;
	printFile_ << "@<TRIPOS>MOLECULE" << endl;
	printFile_ << "name" << endl;
	printFile_ << nAtoms + 1 << "  " << nBonds << endl;
	printFile_ << "SMALL " << endl << "NO_CHARGES" << endl << endl << endl;
	printFile_ << "@<TRIPOS>ATOM" << endl;

	//lantanide
	double zero = 0.0e0;
	printFile_ << " 1 A1 "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< "Au " << endl;

	for (int i = 0; i < nAtoms; i++)
	{
		printFile_ << " " << i + 2 << " A" << i + 1 << "  "
			<< fixed << setprecision(4) << setw(8) << mol0[i].x * scale << "  "
			<< fixed << setprecision(4) << setw(8) << mol0[i].y * scale << "  "
			<< fixed << setprecision(4) << setw(8) << mol0[i].z * scale << "  "
			<< mol0[i].atomlabel << endl;
	}
	printFile_ << "@<TRIPOS>BOND" << endl;
	int k = 1;
	for (size_t i = 0; i < nAtoms; i++)
	{
		printFile_ << " " << k << " 1  " << i + 2 << endl;
		k++;
	}
	for (size_t i = 0; i < bidentateAtomsChosen.size(); i += 2)
	{
		printFile_ << " " << k << " " << bidentateAtomsChosenRotated[i] + 2 << "  "
			<< bidentateAtomsChosenRotated[i + 1] + 2
			<< endl;
		k++;
	}
	printFile_ << "M  END" << endl;
}


void CauchyIndex::identifyIsomer(
	string permutationsFile,
	string filePath,
	string coordinatesFile,
	string countingFile)
{
	IdentifyIsomers identIso_;
	identIso_.coordinatesToPermutation(
		mol0,
		permutationsFile,
		filePath,
		coordinatesFile,
		countingFile);
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

// um arquivo - primeira linha - numero de isomeros - segunda - tipos de atomos - e pronto
// fredmudar --- simasmudar
void CauchyIndex::printAllMoleculesFromFileEnantiomers(string composition)
{
	vector<int> atomTypes;
	vector<int> bidentateChosen;

	// colocar o nome certo aqui.
	readAtomTypesAndBidentateChosenFile(composition + "---atomTypes.txt", atomTypes, bidentateChosen);

	ifstream filePermutations_(("final-" + composition).c_str());
	string line;
	while (!filePermutations_.eof())
	{
		vector<int> permutation = readCauchyNotationsEnantiomers(filePermutations_);
		if (permutation.size() == 0)
			continue;
		string permtString = permutationToString(permutation);
		ofstream fileXyz_((composition + "-" + permtString + ".mol2").c_str());
		printMoleculeMolFormat(permutation, atomTypes, bidentateChosen, fileXyz_);
		fileXyz_.close();
	}
	filePermutations_.close();
#ifdef UNIX
	system(("mkdir " + composition).c_str());
	system(("mv " + composition + "* " + composition).c_str());
#endif

}




// read  new form of isomers.


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

vector<int> CauchyIndex::readCauchyNotationsEnantiomers(ifstream & openendFile_)
{
	int size = mol0.size();
	vector<int> notation;
	if (openendFile_.eof())
		return notation;

	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;
	notation.resize(size);
	stringstream line;
	line << auxline;
	string dummy1, dummy2, dummy3, dummy4;
	line >> dummy1 >> dummy2 >> dummy3 >> dummy4;
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


vector< vector<int> > CauchyIndex::readCauchyNotationsRAMBlockEnantiomers(ifstream & openendFile_, int kInit, int kFinal)
{
	vector< vector<int> > ramBlock;
	vector<int> zero(1);
	zero[0] = -1;
	int k = 1;
	while (true)
	{
		vector<int> isomer = readCauchyNotations(openendFile_);
		if (k >= kInit && k <= kFinal)
		{
			if (isomer.size() == 0)
				ramBlock.push_back(zero);
			else
				ramBlock.push_back(isomer);
		}
		if(k == kFinal)
			break;
		if(isomer.size() != 0)
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

vector<int> CauchyIndex::readNewCauchyNotationsEnantiomers(ifstream & openendFile_)
{
	int size = mol0.size();
	vector<int> notation;
	if (openendFile_.eof())
		return notation;

	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;
	notation.resize(size);

	size_t brack1Temp = auxline.find("]");
	size_t brack1 = auxline.find("[", brack1Temp + 1, 1);
	size_t brack2 = auxline.find("]", brack1Temp + 1, 1);
	string permString = auxline.substr(brack1 + 1, brack2 - brack1 - 1);
	stringstream line;
	line << permString;
	for (size_t i = 0; i < size; i++)
	{
		line >> notation[i];
		notation[i]--;
	}
	return notation;
}


void CauchyIndex::findMissedRotations()
{
	string fileName = "OC-6-Ma2b2c2.csv";

	IsomersToMol isoMol_;
	size_t found = fileName.find("-");
	size_t foundSecond = fileName.find("-", found + 1, 1);
	size_t pointChar = fileName.find(".");
	string composition = fileName.substr(foundSecond + 1, pointChar - foundSecond - 1);
	int nBidentates = 0;
	int coordination = 6;
	ifstream fileIsomers_(fileName.c_str());
	vector<int> atomTypes;
	vector<int> bidentateAtomsChosen;

	isoMol_.readAtomTypesAndBidentateChosenFile(
		fileIsomers_,
		atomTypes,
		bidentateAtomsChosen,
		coordination,
		nBidentates);

	Geometries geomet_;
	double cutAngle;
	vector< vector<int> > allReflections;
	vector<CoordXYZ> molDummy;
	geomet_.geometry6OCReflections(molDummy, cutAngle, allReflections);
	// apply this with bidentates
	// generate S operations by composing reflections with rotations

	ofstream rotations_((fileName + "-rotations.csv").c_str());
	string line;
	rotations_ << "Permutation ; Rotations ";
	while (!fileIsomers_.eof())
	{
		vector<int> permutation = readNewCauchyNotationsEnantiomers(fileIsomers_);
		if (permutation.size() == 0)
		{
			rotations_ << endl;
			continue;
		}
		rotations_ << endl;

		for (size_t i = 0; i < permutation.size(); i++)
			rotations_ << permutation[i] << " ";
		rotations_ << " ; ";

		size_t size = mol0.size();
		vector<int> types1(size);
		vector<int> types2(size);
		for (size_t i = 0; i < size; i++)
		{
			types1[i] = atomTypes[permutation[i]];
			types2[i] = atomTypes[permutation[i]];
		}
		vector<int> bidPos1 = applyPermutationBidentates(permutation, bidentateAtomsChosen);
		if (bidPos1.size() == 0 && bidentateAtomsChosen.size() != 0)
		{
			cout << "findMissedRotations error" << endl;
			exit(1);
		}
		vector<int> bidPos2 = applyPermutationBidentates(permutation, bidentateAtomsChosen);
		if (bidPos2.size() == 0 && bidentateAtomsChosen.size() != 0)
		{
			cout << "findMissedRotations error" << endl;
			exit(1);
		}
		sort(bidPos1.begin(), bidPos1.end());
		sort(bidPos2.begin(), bidPos2.end());
		if (types1 == types2)
		{
			if (bidPos1 == bidPos2)
			{
				rotations_ << " E ; ";
			}
		}
		for (size_t j = 0; j < allRotationTransforms.size(); j++)
		{
			vector<int> auxPerm = applyRotation(permutation, j);
			for (size_t i = 0; i < size; i++)
			{
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
				{
					rotations_ << rotationString(j) << " ; ";
					//cout << "rotation " << rotationString(j) << " preserved" << endl;

					// rotacao preservada
					//return;
				}
			}
		}
	}
	fileIsomers_.close();
}

string CauchyIndex::rotationString(int iRot)
{
	switch (iRot)
	{
	case 0:
		return "C4(0)";
	case 1:
		return "C4-2(0)";
	case 2:
		return "C4-3(0)";
	case 3:
		return "C4(1)";
	case 4:
		return "C4-2(1)";
	case 5:
		return "C4-3(1)";
	case 6:
		return "C4(2)";
	case 7:
		return "C4-2(2)";
	case 8:
		return "C4-3(2)";
	case 9:
		return "C2(0-1)";
	case 10:
		return "C2(0-2)";
	case 11:
		return "C2(0-3)";
	case 12:
		return "C2(0-4)";
	case 13:
		return "C2(1-2)";
	case 14:
		return "C2(1-4)";
	case 15:
		return "C3(0-1-2)";
	case 16:
		return "C3-2(0-1-2)";
	case 17:
		return "C3(0-2-3)";
	case 18:
		return "C3-2(0-2-3)";
	case 19:
		return "C3(0-3-4)";
	case 20:
		return "C3-2(0-3-4)";
	case 21:
		return "C3(0-4-1)";
	case 22:
		return "C3-2(0-4-1)";
	default:
		cout << "rotation on CauchyIndex::rotationString not found" << endl;
		exit(1);
		break;
	}

	return "ERROR";
}















































// 5 - vOC(y)  ;  TBPY (nossa) ; SPY(y)
// 6 - OC (nossa) ; TPR(y)
// 7 - PBPY(y) ; COC (nossa) ; CTPR(y)
// 8 - SAPR (nossa) ; TDD(y - JSD!!!) ; BTPR(y - JBTP!!!) ; HBPY(y) ; CU(y)  
// 9 - CSAPR(y) ; TCTPR (nossa - JTCTPR) ; MFF (cs)
// 10 - JMBIC (nosso) 
// 11 - JCPAPR (nosso)
// 12 - IC (nosso)
void CauchyIndex::setSystem(int system)
{
	vector<double> vectorRotations;
	Geometries geo_;

	vectorRotations = geo_.selectGeometry(system, mol0, cutAngle, reflectionOperation);

	setAllRotations(vectorRotations);
}

