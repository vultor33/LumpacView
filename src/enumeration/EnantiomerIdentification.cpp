#include "EnantiomerIdentification.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "Coordstructs.h"
#include "AuxMath.h"
#include "IsomersToMol.h"


using namespace std;

EnantiomerIdentification::EnantiomerIdentification(){}

EnantiomerIdentification::~EnantiomerIdentification(){}

std::vector<int> EnantiomerIdentification::finalIsomer(
	std::vector<CoordXYZ> &mol0,
	std::vector<int> &atomTypes,
	std::vector<int> &bidentate,
	std::vector<int> &permut1,
	std::vector<int> &permut2,
	std::vector<CoordXYZ> &crystalGeom,
	std::vector<int> &composition,
	std::vector<int> &bidentateCrystal)
{
	vector<int> bidentate1 = bidentate;
	vector<int> bidentate2 = bidentate;
	applyPermutationBidentates(permut1, bidentate1);
	applyPermutationBidentates(permut2, bidentate2);


	//RootMeanSquareDeviation -- oldRmsdToAddressPoints
	vector<int> atomTypes1(atomTypes.size());
	vector<int> atomTypes2(atomTypes.size());
	for (size_t i = 0; i < atomTypes.size(); i++)
	{
		atomTypes1[i] = atomTypes[permut1[i]];
		atomTypes2[i] = atomTypes[permut2[i]];
	}

	int minCount = chooseBestCompareType(composition);
	vector<int> mapping = mapNeededPermutation(minCount, composition, atomTypes1);
	IsomersToMol iso_;
	vector<CoordXYZ> rotGeometry = rotateVectorToZ(mapping[0], crystalGeom);
	double minDist1 = 1.0e99;
	for (size_t i = 1; i < mapping.size(); i++)
	{
		iso_.printCoordXyz(mol0);
		vector<CoordXYZ> rotPermutation = rotateVectorToZ(mapping[i], mol0);
		iso_.printCoordXyz(rotPermutation);
		double auxDist = compareEnantiomersThroughRotation(
			mapping[0],
			mapping[i],
			rotGeometry,
			composition,
			bidentateCrystal,
			rotPermutation,
			atomTypes1,
			bidentate1);
		cout << "mpa1 " << auxDist << endl;
		if (auxDist < minDist1)
			minDist1 = auxDist;
	}
	vector<int> mapping2 = mapNeededPermutation(minCount, composition, atomTypes2);
	double minDist2 = 1.0e99;
	for (size_t i = 1; i < mapping2.size(); i++)
	{
		vector<CoordXYZ> rotPermutation = rotateVectorToZ(mapping2[i], mol0);
		double auxDist = compareEnantiomersThroughRotation(
			mapping2[0],
			mapping2[i],
			rotGeometry,
			composition,
			bidentateCrystal,
			rotPermutation,
			atomTypes2,
			bidentate2);
		cout << "mpa2 " << auxDist << endl;
		if (auxDist < minDist2)
			minDist2 = auxDist;
	}

	if (minDist1 < minDist2)
		return permut1;
	else
		return permut2;
}

std::vector<CoordXYZ> EnantiomerIdentification::rotateVectorToZ(
	int iAtom,
	vector<CoordXYZ> &mol
)
{
	double angle = auxMath_.angleFrom3Points(
		mol[iAtom].x, mol[iAtom].y, mol[iAtom].z,
		0.0e0, 0.0e0, 0.0e0,
		0.0e0, 0.0e0, 1.0e0);

	vector<double> normal = auxMath_.normalVectorFrom3Points(
		mol[iAtom].x, mol[iAtom].y, mol[iAtom].z,
		0.0e0, 0.0e0, 0.0e0,
		0.0e0, 0.0e0, 1.0e0);

	vector< vector<double> > rot = auxMath_.rotationMatrix(normal[0], normal[1], normal[2], angle);

	vector<double> aux;
	vector<CoordXYZ> rotCoord(mol.size());
	for (size_t i = 0; i < mol.size(); i++)
	{
		aux = auxMath_.matrixXVector(rot, mol[i].x, mol[i].y, mol[i].z);
		rotCoord[i].x = aux[0];
		rotCoord[i].y = aux[1];
		rotCoord[i].z = aux[2];
	}

	double angle2 = auxMath_.angleFrom3Points(
		rotCoord[iAtom].x, rotCoord[iAtom].y, rotCoord[iAtom].z,
		0.0e0, 0.0e0, 0.0e0,
		0.0e0, 0.0e0, 1.0e0);

	if (abs(angle2) > 0.01e0)
	{
		cout << "problem on EnantiomerIdentification::rotateVectorToZ"
			<< endl;
		exit(1);

	}
	return rotCoord;
}

int EnantiomerIdentification::chooseBestCompareType(vector<int> &atomTypes1)
{
	std::map<int, unsigned int> rv;
	for (vector<int>::iterator val = atomTypes1.begin(); val != atomTypes1.end(); ++val)
	{
		rv[*val]++;
	}
	int minCount = 1000;
	int iMinCount = -1;
	for (map<int, unsigned int>::iterator count = rv.begin(); count != rv.end(); ++count) 
	{
		if (count->second < minCount)
		{
			iMinCount = count->first;
			minCount = count->second;
		}
	}
	return iMinCount;
}


vector<int> EnantiomerIdentification::mapNeededPermutation(
	int minCount,
	vector<int> &atomTypes1,
	vector<int> &atomTypes2
)
{
	vector<int> permutMap;
	size_t size = atomTypes1.size();
	for (size_t i = 0; i < size; i++)
	{
		if (atomTypes1[i] == minCount)
		{
			permutMap.push_back(i);
			for (size_t j = 0; j < size; j++)
			{
				if (atomTypes2[j] == minCount)
				{
					permutMap.push_back(j);
				}
			}
			break;
		}

	}
	return permutMap;
}


// entrar as duplas e vou querer cuspido o resultado

double EnantiomerIdentification::compareEnantiomersThroughRotation(
	int chosen1,
	int chosen2,
	std::vector<CoordXYZ> &coord1,
	std::vector<int> &atomTypes1,
	std::vector<int> &bidentate1,
	std::vector<CoordXYZ> &coord2,
	std::vector<int> &atomTypes2,
	std::vector<int> &bidentate2)
{
	double minDist = 1.0e99;
	IsomersToMol iso_;//fredmudar
	vector<CoordXYZ> bestCoord;//fredmudar
	vector<int> bestConnect;
	for (int i = 0; i < 63; i++)
	{
		rotateCoordOnChosen(chosen2, coord2, 0.1e0);
		iso_.printCoordXyz(coord2);
		vector< vector<double> > matrDist = generateDistanceMatrix(coord1, coord2);
		vector<int> connect;
		double dist = totalDist(
			chosen1,
			chosen2,
			coord1,
			coord2,
			matrDist,
			atomTypes1,
			bidentate1,
			atomTypes2,
			bidentate2,
			connect);

		if (minDist > dist)
		{
			bestCoord = coord2;//fredmudar
			bestConnect = connect;
			minDist = dist;
		}
	}
	ofstream debug_;
	debug_.open("debug.txt", std::ofstream::out | std::ofstream::app);
	debug_ << "min dist:  " << minDist << endl;
	vector< vector<double> > matrDist = generateDistanceMatrix(coord1, bestCoord);
	debug_ << "connect: ";
	for (size_t i = 0; i < bestConnect.size(); i++)
	{
		debug_ << bestConnect[i] << " ";
	}
	debug_ << endl;

	for (int i = 0; i < matrDist.size(); i++)
	{
		for (int j = 0; j < matrDist.size(); j++)
		{
			debug_ << matrDist[i][j] << "  ";
		}
		debug_ << endl;
	}
	debug_ << endl;
	debug_.close();

	return minDist;
}

double EnantiomerIdentification::totalDist(
	int chosen1,
	int chosen2,
	std::vector<CoordXYZ> &coord1,
	std::vector<CoordXYZ> &coord2,
	std::vector< std::vector<double> > &matrDist,
	std::vector<int> atomTypes1,
	std::vector<int> bidentate1,
	std::vector<int> atomTypes2,
	std::vector<int> bidentate2,
	std::vector<int> &connect)
{
	connect.resize(atomTypes1.size());
	for (size_t i = 0; i < connect.size(); i++)
		connect[i] = -1;
	double totalDist = 0.0e0;

	connect[chosen1] = chosen2;
	totalDist += matrDist[chosen1][chosen2];
	if (find(bidentate1.begin(), bidentate1.end(), chosen1) != bidentate1.end())
	{
		int bidI = searchBidentateMatch(bidentate1, chosen1);
		int bidJ = searchBidentateMatch(bidentate2, chosen2);
		connect[bidI] = bidJ;
		totalDist += matrDist[bidI][bidJ];
		totalDist += auxMath_.norm(
			0.5e0*(coord1[chosen1].x + coord1[bidI].x)
			- 0.5e0*(coord2[chosen2].x + coord2[bidJ].x),
			0.5e0*(coord1[chosen1].y + coord1[bidI].y)
			- 0.5e0*(coord2[chosen2].y + coord2[bidJ].y),
			0.5e0*(coord1[chosen1].z + coord1[bidI].z)
			- 0.5e0*(coord2[chosen2].z + coord2[bidJ].z));
	}

	for (size_t i = 0; i < atomTypes1.size(); i++)
	{
		if (connect[i] != -1)
			continue;

		double minDist = 1.0e99;
		int jMinDist;
		for (size_t j = 0; j < atomTypes1.size(); j++)
		{
			if (find(connect.begin(), connect.end(), j) != connect.end())
				continue;

			if ((minDist > matrDist[i][j]) &&
				(atomTypes1[i] == atomTypes2[j]))
			{
				minDist = matrDist[i][j];
				jMinDist = j;
			}
		}
		connect[i] = jMinDist;
		totalDist += minDist;
		if (find(bidentate1.begin(), bidentate1.end(), i) != bidentate1.end())
		{
			int bidI = searchBidentateMatch(bidentate1, i);
			int bidJ = searchBidentateMatch(bidentate2, jMinDist);
			connect[bidI] = bidJ;
			totalDist += matrDist[bidI][bidJ];
			totalDist += auxMath_.norm(
				0.5e0*(coord1[i].x + coord1[bidI].x)
				- 0.5e0*(coord2[jMinDist].x + coord2[bidJ].x),
				0.5e0*(coord1[i].y + coord1[bidI].y)
				- 0.5e0*(coord2[jMinDist].y + coord2[bidJ].y),
				0.5e0*(coord1[i].z + coord1[bidI].z)
				- 0.5e0*(coord2[jMinDist].z + coord2[bidJ].z));
		}
	}
	return totalDist;
}

vector< vector<double> > EnantiomerIdentification::generateDistanceMatrix(
	vector<CoordXYZ> &coord1,
	vector<CoordXYZ> &coord2)
{
	vector< vector<double> > distMatrix;
	for (size_t i = 0; i < coord1.size(); i++)
	{
		vector<double> auxVec(coord1.size());
		for (size_t j = 0; j < coord1.size(); j++)
		{
			auxVec[j] = auxMath_.norm(
				coord1[i].x - coord2[j].x,
				coord1[i].y - coord2[j].y,
				coord1[i].z - coord2[j].z);
		}
		distMatrix.push_back(auxVec);
	}
	return distMatrix;
}


int EnantiomerIdentification::searchBidentateMatch(
	vector<int> &bidentate,
	int j)
{
	for (size_t i = 0; i < bidentate.size(); i += 2)
	{
		if (bidentate[i] == j)
			return bidentate[i + 1];
		else if (bidentate[i + 1] == j)
			return bidentate[i];
	}
	cout << "On searchBidentateMatch, bidentate no found" << endl;
	exit(1);
	return -1;
}

void EnantiomerIdentification::rotateCoordOnChosen(
	int chosen,
	vector<CoordXYZ> &coord,
	double angle)
{
	vector< vector<double> > rot = auxMath_.rotationMatrix(
		coord[chosen].x, coord[chosen].y, coord[chosen].z,
		angle);

	for (size_t i = 0; i < coord.size(); i++)
	{
		if ((int)i == chosen)
			continue;

		vector<double> auxVec = auxMath_.matrixXVector(rot,
			coord[i].x, coord[i].y, coord[i].z);
		coord[i].x = auxVec[0];
		coord[i].y = auxVec[1];
		coord[i].z = auxVec[2];
	}
}


void EnantiomerIdentification::applyPermutationBidentates(
	const std::vector<int>& permutation,
	std::vector<int>& bidentateAtomsChosen)
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
	bidentateAtomsChosen = bidentateAtomsChosenRotated;
}


