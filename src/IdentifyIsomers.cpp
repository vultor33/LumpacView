#include "IdentifyIsomers.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>

#include "Coordstructs.h"

using namespace std;

IdentifyIsomers::IdentifyIsomers(){}

IdentifyIsomers::~IdentifyIsomers(){}

void IdentifyIsomers::test(vector<CoordXYZ> &mol0)
{
	size_t size = mol0.size();
	vector<int> atomTypes(size);
	atomTypes[0] = 2;
	atomTypes[1] = 0;
	atomTypes[2] = 2;
	atomTypes[3] = 1;
	atomTypes[4] = 0;
	atomTypes[5] = 1;

	string p1 = "0 1 2 3 4 5";
	string file = "test-identify.xyz";
	compareGeometryPermutation(
		atomTypes,
		stringToPermutation(p1,atomTypes.size()),
		mol0,
		file);
	/*
	string p1 = "5 3 1 2 4 6 0";
	string p2;
	p2 = "0 1 2 3 4 5 6";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "0 4 3 2 1 5 6";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "0 2 6 1 4 5 3";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "0 4 1 6 2 5 3";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "3 1 5 4 2 6 0";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "0 1 6 2 5 4 3";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "0 2 6 1 5 4 3";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	p2 = "1 3 2 4 5 0 6";
	compareTwoPermutations(atomTypes, p1, p2, mol0);
	*/

}

void IdentifyIsomers::compareGeometryPermutation(
	vector<int> &atomTypes,
	vector<int> &permutation,
	vector<CoordXYZ> &mol0,
	string fileName)
{
	vector< vector<int> > allT1, allT2;
	vector< vector<double> > allD1, allD2;
	size_t size = atomTypes.size();
	vector<int> atomTypes1(size);
	for (size_t i = 0; i < atomTypes.size(); i++)
	{
		atomTypes1[i] = atomTypes[permutation[i]];
	}
	sortAllDistances(atomTypes1, mol0, allT1, allD1);

	// precisa do tipo e do mol0 MUDAR
	vector<int> composition(atomTypes.size());
	vector<CoordXYZ> geometry = readGeometry(fileName, composition);
	sortAllDistances(composition, geometry, allT2, allD2);

	vector<int> connect;
	vector<double> dist;
	compareIsomers(
		atomTypes1,
		allT1,
		allD1,
		composition,
		allT2,
		allD2,
		connect,
		dist);

	cout << "permutation: ";
	double totalDist = 0.0e0;
	for (size_t i = 0; i < permutation.size(); i++)
	{
		cout << permutation[i] << " ";
		totalDist += dist[i];
	}
	cout << "dist: " << totalDist << endl;


}


void IdentifyIsomers::compareTwoPermutations(
	std::vector<int> &atomTypes,
	std::string &permutation1,
	std::string &permutation2,
	std::vector<CoordXYZ> &mol0)
{
	compareTwoPermutations(
		atomTypes,
		stringToPermutation(permutation1, mol0.size()),
		stringToPermutation(permutation2, mol0.size()),
		mol0);
}


void IdentifyIsomers::compareTwoPermutations(
	vector<int> &atomTypes,
	vector<int> &permutation1,
	vector<int> &permutation2,
	vector<CoordXYZ> &mol0)
{
	vector< vector<int> > allT1, allT2;
	vector< vector<double> > allD1, allD2;
	size_t size = atomTypes.size();
	vector<int> atomTypes1(size);
	vector<int> atomTypes2(size);
	for (size_t i = 0; i < atomTypes.size(); i++)
	{
		atomTypes1[i] = atomTypes[permutation1[i]];
		atomTypes2[i] = atomTypes[permutation2[i]];
	}
	sortAllDistances(atomTypes1, mol0, allT1, allD1);
	sortAllDistances(atomTypes2, mol0, allT2, allD2);
	vector<int> connect;
	vector<double> dist;
	compareIsomers(
		atomTypes1, 
		allT1, 
		allD1, 
		atomTypes2, 
		allT2, 
		allD2,
		connect,
		dist);

	cout << "permutation: ";
	double totalDist = 0.0e0;
	for (size_t i = 0; i < permutation2.size(); i++)
	{
		cout << permutation2[i] << " ";
		totalDist += dist[i];
	}
	cout << "dist: " << totalDist << endl;


}


void IdentifyIsomers::compareIsomers(
	std::vector<int> & atomTypes1,
	std::vector< std::vector<int> > & allSortedTypes1,
	std::vector< std::vector<double> > & allSortedDistances1,
	std::vector<int> & atomTypes2,
	std::vector< std::vector<int> > & allSortedTypes2,
	std::vector< std::vector<double> > & allSortedDistances2,
	std::vector<int> &connections,
	std::vector<double> &conectDifference)
{
	for (size_t i = 0; i < allSortedTypes1.size(); i++)
	{
		vector<int> aT1 = allSortedTypes1[i];
		vector<double> aD1 = allSortedDistances1[i];
		int jMin = -1;
		double difMin = 1.0e99;
		for (size_t j = 0; j < allSortedTypes2.size(); j++)
		{
			if (find(connections.begin(), connections.end(), j) != connections.end())
				continue;

			vector<int> aT2 = allSortedTypes2[j];
			vector<double> aD2 = allSortedDistances2[j];
			double difference = compareTwoAtoms(aT1, aD1, aT2, aD2);
		 	if ((difference < difMin) && (difference > -0.1e0))
			{
				jMin = j; // se ja tem algum q foi minimo morre aqui.
				difMin = difference;
			}
		}
		connections.push_back(jMin);
		conectDifference.push_back(difMin);
	}
}


double IdentifyIsomers::compareTwoAtoms(
	std::vector<int> aT1,
	std::vector<double> aD1,
	std::vector<int> aT2,
	std::vector<double> aD2)
{
	//ver se � o mesmo tipo.
	double difference = 0.0e0;
	vector<int> test1 = aT1;
	vector<int> test2 = aT2;
	sort(test1.begin(), test1.end());
	sort(test2.begin(), test2.end());
	if (test1 != test2)
		return -1.0e0;
	else
	{
		while (aT1.size() != 0)
		{
			int kShift = 0;
			for (size_t j = 0; j < aT1.size(); j++)
			{
				if ((j + kShift) > -1 + aT2.size())
					kShift = 0;
				if (aT1[j] == aT2[j + kShift])
				{
					difference += abs(aD1[j] - aD2[j + kShift]);
					aT1.erase(aT1.begin() + j);
					aD1.erase(aD1.begin() + j);
					aT2.erase(aT2.begin() + j + kShift);
					aD2.erase(aD2.begin() + j + kShift);
					j--;
				}
				else
				{
					for (size_t k = j; k < aT2.size(); k++)
					{
						if (aT1[j] == aT2[j + k])
						{
							kShift = k;
							j--;
							break;
						}
						if (k == aT2.size() - 1)
						{
							kShift = 0;
						}
					}
				}
			}
		}
	}
	return difference;
}


void IdentifyIsomers::sortAllDistances(
	std::vector<int> &atomTypes,
	std::vector<CoordXYZ> &mol0,
	std::vector< std::vector<int> > &allSortedTypes,
	std::vector< std::vector<double> > &allSortedDistances)
{
	for (int k = 0; k < mol0.size(); k++)
	{
		vector<int> tempTypes;
		vector<double> tempDistance;
		sortDistancesK(k, atomTypes, mol0, tempTypes, tempDistance);
		allSortedTypes.push_back(tempTypes);
		allSortedDistances.push_back(tempDistance);
	}


}

void IdentifyIsomers::calculateDistancesK(
	int k,
	vector<CoordXYZ> &mol0,
	vector<double> &distances)
{
	double auxDist;
	for (size_t i = 0; i < mol0.size(); i++)
	{
		if (i == k)
			continue;
		auxDist = sqrt(
			(mol0[k].x - mol0[i].x)*(mol0[k].x - mol0[i].x)
			+ (mol0[k].y - mol0[i].y)*(mol0[k].y - mol0[i].y)
			+ (mol0[k].z - mol0[i].z)*(mol0[k].z - mol0[i].z));
		distances.push_back(auxDist);
	}
}

void IdentifyIsomers::sortDistancesK(
	int k,
	vector<int> &atomTypes,
	vector<CoordXYZ> &mol0,
	vector<int> &sortedTypes,
	vector<double> &sortedDistances)
{
	// ESQUECA BIDENTADO NESTE MOMENTO
	vector<double> tempDistances;
	//se os grupos nao baterem o que eu faco?
	calculateDistancesK(k, mol0, tempDistances);
	vector<int> reducedTypes = atomTypes;
	reducedTypes.erase(reducedTypes.begin() + k);
	while (reducedTypes.size() > 0)
	{
		int smallest = indexofSmallestElement(tempDistances);
		//add to new vectors
		int smallestType = reducedTypes[smallest];
		sortedDistances.push_back(tempDistances[smallest]);
		sortedTypes.push_back(reducedTypes[smallest]);
		//erase
		tempDistances.erase(tempDistances.begin() + smallest);
		reducedTypes.erase(reducedTypes.begin() + smallest);
		insertSortedType(smallestType, reducedTypes, tempDistances, sortedTypes, sortedDistances);
	}

}

void IdentifyIsomers::insertSortedType(
	int type,
	vector<int> & reducedTypes,
	vector<double> & distances,
	vector<int> & newTypes,
	vector<double> & newDistances)
{
	vector<double> thisTypeDistances;
	vector<int> thisTypes;
	for (size_t i = 0; i < reducedTypes.size(); i++)
	{
		if (reducedTypes[i] == type)
		{
			thisTypeDistances.push_back(distances[i]);
			thisTypes.push_back(reducedTypes[i]);
			reducedTypes.erase(reducedTypes.begin() + i);
			distances.erase(distances.begin() + i);
			i--;
		}
	}

	// sort thisTypeDistances and add it.
	if (thisTypeDistances.size() > 1)
		sort(thisTypeDistances.begin(), thisTypeDistances.end());

	if (thisTypeDistances.size() != 0)
	{
		newDistances.insert(newDistances.end(), thisTypeDistances.begin(), thisTypeDistances.end());
		newTypes.insert(newTypes.end(), thisTypes.begin(), thisTypes.end());
	}
}


int IdentifyIsomers::indexofSmallestElement(vector<double> & vector)
{
	int index = 0;
	int size = vector.size();
	for (int i = 1; i < size; i++)
	{
		if (vector[i] < vector[index])
			index = i;
	}
	return index;
}


vector<int> IdentifyIsomers::stringToPermutation(string entryString, size_t size)
{
	stringstream convert;
	convert << entryString;
	vector<int> permutation(size);
	for (size_t i = 0; i < size; i++)
	{
		convert >> permutation[i];
	}
	return permutation;




}

std::vector<CoordXYZ> IdentifyIsomers::readGeometry(std::string fileName, std::vector<int> &composition)
{
	ifstream readG_(fileName.c_str());
	string line;
	getline(readG_, line);
	getline(readG_, line);
	vector<CoordXYZ> geometry(composition.size());
	for(size_t i = 0; i < composition.size(); i++)
	{
		getline(readG_, line);
		stringstream convert;
		convert << line;
		convert >> composition[i]
			>> geometry[i].x
			>> geometry[i].y
			>> geometry[i].z;

		double r = sqrt(
			(geometry[i].x * geometry[i].x)
			+ (geometry[i].y * geometry[i].y)
			+ (geometry[i].z * geometry[i].z));
	}
	readG_.close();
	return geometry;



}