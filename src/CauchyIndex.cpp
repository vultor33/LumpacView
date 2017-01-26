#include "CauchyIndex.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <time.h>

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

	cout << "estou no fim" << endl;

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
	string cauchyFileName = "allCauchyRotatations.txt";
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
			
			cout << "i:  " << i << " j: " << j << "  dist: " << atomPosition << endl;
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
			<< mol[i].x << "  "
			<< mol[i].y << "  "
			<< mol[i].z << endl;
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

std::vector<int> CauchyIndex::applyPermutation(
	const std::vector<int> & permutation, 
	const std::vector<std::string> & atoms,
	const std::vector<int> & bidentateAtomsChosen)
{
	size_t size = atoms.size();
	vector<int> bidentateAtomsChosenRotated = bidentateAtomsChosen;
	for (size_t i = 0; i < size; i++)
	{
		mol0[i].atomlabel = atoms[permutation[i]];
	}
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

	vector<int> bidentateAtomsChosenRotated = applyPermutation(permutation, atoms, bidentateAtomsChosen);

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
	for (int i = 0; i < nAtoms; i++)
	{
		if (i < (int)mol0.size())
		{
			printFile_ << mol0[i].atomlabel << "  "
				<< mol0[i].x << "  "
				<< mol0[i].y << "  "
				<< mol0[i].z << endl;
		}
		else
		{
			printFile_ << bidentates[i - mol0.size()].atomlabel << "  "
				<< bidentates[i - mol0.size()].x << "  "
				<< bidentates[i - mol0.size()].y << "  "
				<< bidentates[i - mol0.size()].z << endl;
		}
	}
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
	vector<int> notation;
	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;
	notation.resize(systemSize);
	stringstream line;
	line << auxline;
	for (size_t i = 0; i < systemSize; i++)
	{
		line >> notation[i];
	}
	return notation;
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
		cutAngle = auxMath_._pi / 2.0e0;
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
		mol0[1].y = 0.000;
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
		// add rotations
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
