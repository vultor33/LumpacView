#ifndef FINDISOMERS_H
#define FINDISOMERS_H

#include <fstream>
#include <sys/stat.h>
//#include <unistd.h>
#include <string>

#include "Coordstructs.h"
#include "Ligand.h"

class FindIsomers
{
public:
	FindIsomers();

	~FindIsomers();

	void start();

private:
	int counter;

	std::vector< std::string > inputInformations;

	double identicalStructuresLimit;

	void permutation(int nMax);

	void ligandFilePositionPermutation(std::vector<int> & permutation);

	unsigned int factorial(unsigned int n);

	std::vector< Ligand > setThisPermutationLig(std::vector<int> permutation, std::vector<Ligand> & ligOriginal);

	std::vector< CoordXYZ > setThisPermutationAtoms(std::vector<int> permutation, std::vector<CoordXYZ> &  originAtoms);

	std::vector<CoordXYZ> ligandToCoordXYZ(std::vector<Ligand> & allLigands);

	std::vector<CoordXYZ> readMidXyz(std::ifstream & openStream);

	bool doOverlayWithPreviousConfigurations(std::vector<CoordXYZ> & atomsPointPermutation);

	std::string permutationToString(std::vector<int> & permutation);

	void appendPrintCoordXYZ(std::vector<CoordXYZ> & allAtoms, std::string fName, std::string title);

	void appendPrintCoordXYZ(std::vector<Ligand> & allAtoms, std::string fName, std::string title);

	inline bool exists_test0(const std::string& name)
	{
		std::ifstream f(name.c_str());
		return f.good();
	}

	std::string fileAllIsomers;

	std::ofstream streamAllIsomers_;



};

#endif
