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
	void findMapToReferencePermutation(int filePermutation, std::vector< std::vector<int> > & allPerm, std::vector< std::vector<int> > &internalPerm, int & mapToReferenceI, double & mapToReferenceRms);

	void printSupersition(int lowestFilePosition, int lowestPermutationPosition, std::vector< std::vector<int> > & allPerm, std::vector< std::vector<int> > & internalPerm);

	std::vector< std::string > setThisPermutation(std::vector<int> permutation);
	
	std::vector< Ligand > setThisPermutationLig(std::vector<int> permutation, std::vector<Ligand> & ligOriginal);

	void printCoordXYZ(std::vector<CoordXYZ> & allAtoms, std::string fName, std::string title = "useless");

	std::vector<CoordXYZ> ligandToCoordXYZ(std::vector<Ligand> & allLigands);

	unsigned int factorial(unsigned int n);
	
	//nova estrutura
	void ligandPointPositionPermutaion(int nMax);
	void ligandFilePositionPermutation(std::vector<int> & permutation);
	std::vector< Ligand > bestLigand;
	double bestRmsd;
	/////////////////



	std::vector< std::vector<int> > allFactorialPermutations(const int nMax);

	std::vector< std::string > originalPermutation;

	std::string referenceFile;
};

#endif