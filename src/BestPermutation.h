#ifndef BESTPERMUTATION_H
#define BESTPERMUTATION_H

#include <vector>
#include <string>

#include "Coordstructs.h"
#include "Ligand.h"

class BestPermutation
{
public:
	BestPermutation(std::vector< std::string > originalPermutation_in, std::string referenceFile);

	~BestPermutation();

	void findBestPermutation();

private:
	std::vector< std::string > setThisPermutation(std::vector<int> permutation);

	void printCoordXYZ(std::vector<CoordXYZ> & allAtoms, std::string fName);

	std::vector<CoordXYZ> ligandToCoordXYZ(std::vector<Ligand> & allLigands);

	unsigned int factorial(unsigned int n);

	std::vector< std::vector<int> > allFactorialPermutations(const int nMax);

	std::vector< std::string > originalPermutation;

	std::string referenceFile;
};

#endif