#ifndef BESTPERMUTATION_H
#define BESTPERMUTATION_H

#include <vector>
#include <string>

#include "Coordstructs.h"
#include "Ligand.h"

class BestPermutation
{
public:
	BestPermutation(std::string referenceFile);

	~BestPermutation();

	void findBestPermutation();

private:
	//data
	std::vector< Ligand > bestLigand;
	double bestRmsd;
	std::string referenceFile;
	
	std::vector< Ligand > setThisPermutationLig(std::vector<int> permutation, std::vector<Ligand> & ligOriginal);

	void printCoordXYZ(std::vector<CoordXYZ> & allAtoms, std::string fName, std::string title = "useless");

	std::vector<CoordXYZ> ligandToCoordXYZ(std::vector<Ligand> & allLigands);
	
	void ligandPointPositionPermutaion(int nMax);

	void ligandFilePositionPermutation(std::vector<int> & permutation);

	std::vector< std::vector<int> > allFactorialPermutations(const int nMax);

	unsigned int factorial(unsigned int n);
};

#endif