#ifndef IDENTIFYISOMERS_H
#define IDENTIFYISOMERS_H

#include <vector>
#include "Coordstructs.h"

class IdentifyIsomers
{
public:
	IdentifyIsomers();
	~IdentifyIsomers();

	bool compareIsomers(
		std::vector<int> & atomTypes1,
		std::vector< std::vector<int> > & allSortedTypes1,
		std::vector< std::vector<double> > & allSortedDistances1,
		std::vector<int> & atomTypes2,
		std::vector< std::vector<int> > & allSortedTypes2,
		std::vector< std::vector<double> > & allSortedDistances2);


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

private:
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

};

#endif
