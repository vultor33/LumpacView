#include "IdentifyIsomers.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "Coordstructs.h"
#include "IsomersToMol.h"

using namespace std;

IdentifyIsomers::IdentifyIsomers(){}

IdentifyIsomers::~IdentifyIsomers(){}

void IdentifyIsomers::coordinatesToPermutation(
	vector<CoordXYZ> &mol0,
	string permutationsFile,
	string coordinatesFile)
{
	IsomersToMol isoMol_;
	vector<int> atomTypes;
	vector<int> bidentateChosen;
	vector<string> allPerm = isoMol_.readAllPermutations(
		permutationsFile, 
		atomTypes,
		bidentateChosen);

	size_t size = mol0.size();
	string file = coordinatesFile;
	vector<int> composition(atomTypes.size());
	vector<int> bidentateGeometry;
	vector<CoordXYZ> outGeometry = readGeometry(
		file, 
		composition, 
		bidentateChosen.size(),
		bidentateGeometry);

	double minimumDist = 1.0e99;
	int minimumPermut = -1;
	double auxDist;
	for (size_t i = 0; i < allPerm.size(); i++)
	{
		double auxDist = compareGeometryPermutation(
			atomTypes,
			stringToPermutation(allPerm[i], atomTypes.size()),
			mol0,
			bidentateChosen,
			composition,
			outGeometry,
			bidentateGeometry);
		if (auxDist < minimumDist)
		{
			minimumDist = auxDist;
			minimumPermut = i;
		}
	}

	isoMol_.printSingleMol(
		stringToPermutation(allPerm[minimumPermut], atomTypes.size()),
		atomTypes,
		bidentateChosen,
		coordinatesFile);

	isoMol_.printSingleMol(
		stringToPermutation("4 3 2 8 1 5 7 0 6",atomTypes.size()),
		atomTypes,
		bidentateChosen,
		coordinatesFile);

	
}

double IdentifyIsomers::compareGeometryPermutation(
	vector<int> &atomTypes,
	vector<int> &permutation,
	vector<CoordXYZ> &mol0,
	vector<int> &bidentatePermutation,
	vector<int> &composition,
	vector<CoordXYZ> &outGeometry,
	vector<int> &bidentateGeometry)
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

	sortAllDistances(composition, outGeometry, allT2, allD2);

	vector<int> connect;
	vector<int> bidentatePermutationRotated = bidentatePermutation;
	applyPermutationBidentates(permutation, bidentatePermutationRotated);
	vector<double> dist;
	compareIsomers(
		atomTypes1,
		bidentatePermutationRotated,
		allT1,
		allD1,
		composition,
		bidentateGeometry,
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
	cout << "dist: " << setprecision(16) << totalDist << endl;
	return totalDist;


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
	cout << "not working, need bidentates" << endl;
	exit(1);
	/*
	compareIsomers(
		atomTypes1, 
		allT1, 
		allD1, 
		atomTypes2, 
		allT2, 
		allD2,
		connect,
		dist);
	*/
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
	std::vector<int> & bidentates1,
	std::vector< std::vector<int> > & allSortedTypes1,
	std::vector< std::vector<double> > & allSortedDistances1,
	std::vector<int> & atomTypes2,
	std::vector<int> & bidentates2,
	std::vector< std::vector<int> > & allSortedTypes2,
	std::vector< std::vector<double> > & allSortedDistances2,
	std::vector<int> &connections,
	std::vector<double> &conectDifference)
{
	connections.resize(atomTypes1.size()); //fredmudar
	conectDifference.resize(atomTypes1.size());
	for (size_t i = 0; i < connections.size(); i++)
		connections[i] = -1;
	for (size_t i = 0; i < allSortedTypes1.size(); i++)
	{
		if (connections[i] != -1)
			continue;

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
		// a pergunta e - o i tem um bidentado? se sim, qual e?
		// i = 3;
		// o j tambem tem um bidentado.
		connections[i] = jMin;
		conectDifference[i] = difMin;
		if (find(bidentates1.begin(), bidentates1.end(), i) != bidentates1.end())
		{
			int bidI = searchBidentateMatch(bidentates1, i);
			int bidJ = searchBidentateMatch(bidentates2, jMin);
			connections[bidI] = bidJ;
			vector<int> aT2n = allSortedTypes2[bidJ];
			vector<double> aD2n = allSortedDistances2[bidJ];
			vector<int> aT1n = allSortedTypes1[bidI];
			vector<double> aD1n = allSortedDistances1[bidI];
			conectDifference[bidI] = compareTwoAtoms(aT1n, aD1n, aT2n, aD2n);	
		}
	}
}


double IdentifyIsomers::compareTwoAtoms(
	std::vector<int> aT1,
	std::vector<double> aD1,
	std::vector<int> aT2,
	std::vector<double> aD2)
{
	//ver se é o mesmo tipo.
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

std::vector<CoordXYZ> IdentifyIsomers::readGeometry(
	std::string fileName, 
	std::vector<int> &composition,
	int bidSize,
	std::vector<int> &bidentateGeometry)
{
	ifstream readG_(fileName.c_str());
	string line;
	getline(readG_, line);
	getline(readG_, line);
	vector<CoordXYZ> geometry(composition.size());
	for (size_t i = 0; i < composition.size(); i++)
	{
		getline(readG_, line);
		stringstream convert;
		convert << line;
		string comp;
		convert >> comp
			>> geometry[i].x
			>> geometry[i].y
			>> geometry[i].z;
		composition[i] = reverseComposition(comp);
	}
	for (int i = 0; i < bidSize; i++)
	{
		getline(readG_, line);
		if (line == "")
		{
			cout << "Geometry bidentates not found" << endl;
			exit(1);
		}
		stringstream convert;
		convert << line;
		int auxBid;
		convert >> auxBid;
		bidentateGeometry.push_back(auxBid);
	}

	translateToCenterOfMassAndReescale(geometry);

	readG_.close();
	return geometry;
}

int IdentifyIsomers::reverseComposition(string comp)
{
	if (comp == "Ca")
		return 0;
	else if (comp == "O")
		return 1;
	else if (comp == "H")
		return 2;
	else if (comp == "Mg")
		return 3;
	else if (comp == "Be")
		return 4;
	else if (comp == "I")
		return 5;
	else if (comp == "Cl")
		return 6;
	else if (comp == "Br")
		return 7;
	else if (comp == "Na")
		return 8;
	else if (comp == "He")
		return 9;
	else if (comp == "N")
		return 10;
	else if (comp == "C")
		return 11;

	return 100;
}

void IdentifyIsomers::translateToCenterOfMassAndReescale(vector<CoordXYZ> &coord)
{
	double cmx = 0.0e0;
	double cmy = 0.0e0;
	double cmz = 0.0e0;
	for (size_t i = 0; i < coord.size(); i++)
	{
		cmx += coord[i].x;
		cmy += coord[i].y;
		cmz += coord[i].z;
	}
	cmx /= (double)coord.size();
	cmy /= (double)coord.size();
	cmz /= (double)coord.size();

	for (size_t i = 0; i < coord.size(); i++)
	{
		coord[i].x -= cmx;
		coord[i].y -= cmy;
		coord[i].z -= cmz;
		double r = sqrt(
			(coord[i].x * coord[i].x)
			+ (coord[i].y * coord[i].y)
			+ (coord[i].z * coord[i].z));
		coord[i].x /= r;
		coord[i].y /= r;
		coord[i].z /= r;
	}
}

bool IdentifyIsomers::compareConnectBidentates(
	std::vector<int> &connect,
	std::vector<int> &bidPermut,
	std::vector<int> &bidGeom)
{
	for (size_t i = 0; i < bidPermut.size(); i+=2)
	{
		int c1, c2;
		c1 = connect[bidPermut[i]];
		c2 = connect[bidPermut[i + 1]];

		//c1 c2 are bidentates?
		bool areBidentates = false; 
		for (size_t j = 0; j < bidGeom.size(); j += 2)
		{
			if ((c1 == bidGeom[j]) && (c2 == bidGeom[j + 1]))
			{
				areBidentates = true;
				break;
			}
			else if ((c2 == bidGeom[j]) && (c1 == bidGeom[j + 1]))
			{
				areBidentates = true;
				break;
			}
		}
		if (!areBidentates)
			return false;
	}
	return true;
}

int IdentifyIsomers::searchBidentateMatch(
	vector<int> &bidentate,
	int j)
{
	for (size_t i = 0; i < bidentate.size(); i+= 2)
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

void IdentifyIsomers::applyPermutationBidentates(
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


/*
string p1;
p1 = "0 1 2 3 4 5";
compareGeometryPermutation(
atomTypes,
stringToPermutation(p1,atomTypes.size()),
mol0,
composition,
outGeometry);
p1 = "0 1 4 3 2 5";
compareGeometryPermutation(
atomTypes,
stringToPermutation(p1, atomTypes.size()),
mol0,
composition,
outGeometry);
p1 = "0 1 3 4 2 5";
compareGeometryPermutation(
atomTypes,
stringToPermutation(p1, atomTypes.size()),
mol0,
composition,
outGeometry);
p1 = "0 1 3 5 4 2";
compareGeometryPermutation(
atomTypes,
stringToPermutation(p1, atomTypes.size()),
mol0,
composition,
outGeometry);
p1 = "0 2 3 4 5 1";
compareGeometryPermutation(
atomTypes,
stringToPermutation(p1, atomTypes.size()),
mol0,
composition,
outGeometry);
p1 = "0 1 2 4 3 5";
compareGeometryPermutation(
atomTypes,
stringToPermutation(p1, atomTypes.size()),
mol0,
composition,
outGeometry);
*/







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
