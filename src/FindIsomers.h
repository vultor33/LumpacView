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

	std::vector<CoordXYZ> generateSelectedIsomer(std::vector<int> permutation, std::vector<std::string > inputInformations);

	void printSelectedIsomer(std::vector<int> permutation, std::vector<std::string > inputInformations, std::string outputName);

	std::vector<CoordXYZ> buildComplexWithSelectedIsomer(
		std::vector<int> & permutation, 
		std::string projectName, 
		std::string method, 
		std::vector< std::string > & options, 
		std::string & mopacExecPath);

	void buildComplexWithSelectedIsomer();

private:
	//data
	int counter;
	int allCounter;
	int maxCounter;
	double wall0;
	double wall1;
	std::vector<int> bidentateAtoms;
	std::vector<int> bidentateAtomsCombination;
	double bidentateAngleCut;

	std::vector< std::vector<CoordXYZ> > allConfigurations;

	std::vector< std::vector<int> >  allConfigurationPermutation;

	bool useFile;

	std::vector< std::string > inputInformations;

	double identicalStructuresLimit;

	std::string fileAllIsomers;

	std::ofstream streamAllIsomers_;

	//functions
	void permutation(int nMax);

	void ligandFilePositionPermutation(std::vector<int> & permutation);

	unsigned int factorial(unsigned int n);

	std::vector< Ligand > setThisPermutationLig(std::vector<int> permutation, std::vector<Ligand> & ligOriginal);

	std::vector< CoordXYZ > setThisPermutationAtoms(std::vector<int> permutation, std::vector<CoordXYZ> &  originAtoms);

	void aplyPermutationBidentate(std::vector<int> permutation, std::vector<CoordXYZ> & atomsPointPermutation);

	std::vector<CoordXYZ> readMidXyz(std::ifstream & openStream);

	bool doOverlayWithPreviousConfigurations(std::vector<CoordXYZ> & atomsPointPermutation);

	std::string permutationToString(std::vector<int> & permutation);

	void appendPrintCoordXYZ(std::vector<CoordXYZ> & allAtoms, std::string fName, std::string title);

	void appendPrintCoordXYZ(std::vector<CoordXYZ> & allAtoms, std::ofstream & fName_, std::string title);

	void appendPrintCoordXYZ(std::vector<Ligand> & allAtoms, std::string fName, std::string title);

	inline bool exists_test0(const std::string& name)
	{
		std::ifstream f(name.c_str());
		return f.good();
	}

	std::vector<CoordXYZ> ligandToCoordXYZ(std::vector<Ligand> & allLigands);

	void readInput();

	void readInputpermutation(
		std::vector< std::string > & inputPermutations, 
		std::vector<int> & permutation,
		std::vector< std::string > & options,
		std::string & execpath
		);

};

#endif


/*
VOU PRECISAR DE UM FLAG FALANDO QUAL COMBINACAO E O BIDENTADO





*/