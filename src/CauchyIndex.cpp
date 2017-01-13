#include "CauchyIndex.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "RootMeanSquareDeviation.h"
#include "Coordstructs.h"
#include "AuxMath.h"

using namespace std;

CauchyIndex::CauchyIndex()
{
	bidentateLabels.resize(6);
	bidentateLabels[0] = "He";
	bidentateLabels[1] = "Ne";
	bidentateLabels[2] = "Ar";
	bidentateLabels[3] = "Kr";
	bidentateLabels[4] = "Xe";
	bidentateLabels[5] = "Rn";
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


vector<int> CauchyIndex::getCauchy(string fName) 
{
//	RootMeanSquareDeviation rmsd_;
	return vector<int>();
}

void CauchyIndex::calculateAllIndexes()
{
	setSystem(5);

	vector< vector<int> > allIndexes;

	for (size_t i = 0; i < allRotations.size(); i++)
	{
		vector<int> cauchyI = getCauchy(i);

		allIndexes.push_back(cauchyI);
	}

	calculateBidentateMap();
}


void CauchyIndex::getAllRotationMatrix()
{
	setSystem(5);
	vector< vector<int> > allIndexRotations;
	for (size_t i = 0; i < allRotations.size(); i++)
	{
		vector<int> rotI = getCauchy(i);
		allIndexRotations.push_back(rotI);
	}
}


std::vector<int> CauchyIndex::getCauchy(int rotation)
{
	setSystem(5);
	
	vector<CoordXYZ> mol = mol0;

	vector< vector<double> > rot = allRotations[rotation].mRot;

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
	allRotations.resize(allRotationsVector.size() / 4);
	for (size_t i = 0; i < allRotationsVector.size(); i += 4)
	{
		allRotations[i / 4].mRot = auxMath_.rotationMatrix(
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
	double cutAngle = 3.0e0 * auxMath_._pi / 4.0e0;
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

void CauchyIndex::printMolecule(
	vector<int> & permutation,
	vector<string> & atoms,
	vector<int> bidentateAtomsChosen)
{
	for (size_t i = 0; i < atoms.size(); i++)
		mol0[i].atomlabel = atoms[permutation[i]];

	vector<CoordXYZ> bidentates;
	int k = 0;
	for (size_t i = 0; i < bidentateAtomsChosen.size(); i+=2)
	{
		int mapPosition = bidentateMap[i][i + 1];
		if (mapPosition == -1)
			cout << "CauchyIndex::printMolecule - wrong bidentate choice" << endl;

		CoordXYZ auxBidentateAtom = molBidentate[mapPosition];
		auxBidentateAtom.atomlabel = bidentateLabels[k];
		k++;
		bidentates.push_back(auxBidentateAtom);
	}

	int nAtoms = mol0.size() + bidentates.size();
	ofstream of_("printCauchyMolecule.xyz");
	of_ << nAtoms << endl << "title" << endl;
	for (int i = 0; i < nAtoms; i++)
	{
		if (i < mol0.size())
		{
			of_ << mol0[i].atomlabel << "  "
				<< mol0[i].x << "  "
				<< mol0[i].y << "  "
				<< mol0[i].z << endl;
		}
		else
		{
			of_ << bidentates[i - mol0.size()].atomlabel << "  "
				<< bidentates[i - mol0.size()].x << "  "
				<< bidentates[i - mol0.size()].y << "  "
				<< bidentates[i - mol0.size()].z << endl;
		}
	}

	of_.close();


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




