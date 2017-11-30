#ifndef IDENTIFYISOMERS_H
#define IDENTIFYISOMERS_H

#include <vector>
#include <string>
#include "Coordstructs.h"

class IdentifyIsomers
{
public:
	IdentifyIsomers();
	~IdentifyIsomers();

	void test(std::vector<CoordXYZ> & mol0);

	void compareGeometryPermutation(
		std::vector<int> &atomTypes,
		std::vector<int> &permutation,
		std::vector<CoordXYZ> &mol0,
		std::vector<int> &composition,
		std::vector<CoordXYZ> & outGeometry);

	void compareTwoPermutations(
		std::vector<int> &atomTypes,
		std::string &permutation1,
		std::string &permutation2,
		std::vector<CoordXYZ> &mol0);

private:
	void compareTwoPermutations(
		std::vector<int> &atomTypes,
		std::vector<int> &permutation1,
		std::vector<int> &permutation2,
		std::vector<CoordXYZ> &mol0);

	int indexofSmallestElement(std::vector<double> &vector);

	void insertSortedType(
		int type,
		std::vector<int> & reducedTypes,
		std::vector<double> & distances,
		std::vector<int> & newTypes,
		std::vector<double> & newDistances);

	double compareTwoAtoms(
		std::vector<int> aT1,
		std::vector<double> aD1,
		std::vector<int> aT2,
		std::vector<double> aD2);


	void compareIsomers(
		std::vector<int> & atomTypes1,
		std::vector< std::vector<int> > & allSortedTypes1,
		std::vector< std::vector<double> > & allSortedDistances1,
		std::vector<int> & atomTypes2,
		std::vector< std::vector<int> > & allSortedTypes2,
		std::vector< std::vector<double> > & allSortedDistances2,
		std::vector<int> &connections,
		std::vector<double> &difference);


	void sortAllDistances(
		std::vector<int> &atomTypes,
		std::vector<CoordXYZ> &mol0,
		std::vector< std::vector<int> > &allSortedTypes,
		std::vector< std::vector<double> > &allSortedDistances);

	void calculateDistancesK(
		int k,
		std::vector<CoordXYZ> &mol0,
		std::vector<double> &distances);

	void sortDistancesK(
		int k,
		std::vector<int> &atomTypes,
		std::vector<CoordXYZ> &mol0,
		std::vector<int> &sortedTypes,
		std::vector<double> &sortedDistances);

	std::vector<int> stringToPermutation(std::string entryString, size_t size);

	std::vector<CoordXYZ> readGeometry(std::string fileName, std::vector<int> &composition);

};

#endif
