#include "CauchyIndex.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

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
	int nMax = mol0.size();
	long int size = factorial(nMax); // 10+ explosion

	size_t nRotations = allRotationTransforms.size() + 1;
	vector<liPermutation> allPermutations;

	int * myints;
	myints = new int[nMax];
	for (int i = 0; i < nMax; i++)
		myints[i] = i;
	std::sort(myints, myints + nMax);
	vector<int> permutation(nMax);
	bool equal;
	do
	{
		equal = false;
		for (int i = 0; i < nMax; i++)
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

	} while (std::next_permutation(myints, myints + nMax));

	cout << "estou no fim" << endl;

	delete[] myints;
}



void CauchyIndex::rotationTest(
	vector<string> & atoms, 
	vector<int> & bidentateAtomsChosen)
{
	size_t size = atoms.size();
	vector<int> permutation(size);
	for (size_t i = 0; i < size; i++)
		permutation[i] = i;

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

void CauchyIndex::setSystem(int system)
{
	vector<double> vectorRotations;

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
		vectorRotations.resize(mol0.size() * 4);
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

	default:
		cout << "CauchyIndex::setSystem - system not found" << endl;
		exit(1);
		break;
	}
	setAllRotations(vectorRotations);
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




