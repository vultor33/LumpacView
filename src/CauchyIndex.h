#ifndef CAUCHYINDEX_H
#define CAUCHYINDEX_H

#include <vector>
#include <string>
#include <fstream>

#include "Coordstructs.h"
#include "AuxMath.h"

/*
WARNING - new systems need to be set on "setSystem"
          need symmetry rotations and, for bidentate: cutAngle.
*/


struct cauchyRotation
{
	std::vector< std::vector<double> > mRot;
};

struct liPermutation
{
	std::vector< std::vector<int> > rotPermutations;

};

class CauchyIndex
{
public:
	CauchyIndex(int iSystem);
	
	~CauchyIndex();

	//print with mol0 and molBidentate
	// COLOCAR CONST EM TUDO
	void printMolecule(
		std::vector<int> & permutation,
		const std::vector<std::string> & atoms,
		const std::vector<int> & bidentateAtomsChosen,
		std::ofstream & printFile_);

	void rotationTest(
		std::vector<std::string> &atoms, 
		std::vector<int> & bidentateAtomsChosen);

	void generateAllIndependentIsomers();

	void generateAllIndependentIsomersIO();

	void generateAllIndependentIsomersRuntimeRotations();

	void generateAllIndependentIsomersRuntimeRotationsAndReadBlock(std::string blockFileName);

	void generateAllIndependentIsomersWithFlag(std::string blockFileName, std::string code);

	void generateAllIndependentIsomers12(std::string blockFileName);

	//0-true ; 1-false ; 2-bidentate problem
	int compareTwoIsomersWithLabels(
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateAtomsChosen,
		std::vector<int> & permutationIsomer1,
		std::vector<int> & permutationIsomer2
		);

	// true = equal || cant use the same permutation
	bool compareTwoIsomers(
		std::vector<int> & permutationIsomer1,
		std::vector<int> & permutationIsomer2
		);

	void printBlock(int nPieces);

	void mergeBlocks(std::vector<std::string> & allBlockNames, int nPieces);

	std::vector<CoordXYZ> getPoints();

private:
	//data rotations
	std::vector<CoordXYZ> mol0;
	std::vector<CoordXYZ> molBidentate; // all bidentate possible positions marked with atoms.
	std::vector< std::vector<int> > bidentateMap; // between i and j we have k. K position is the same shown on molBidentate.
	std::vector<cauchyRotation> allRotationsMatrix;
	std::vector< std::vector<int> > allRotationTransforms;
	std::vector<std::string> bidentateLabels;
	double cutAngle;

	//data all isomers

	// set system
	void calculateAllIndexes(int iSystem);
	void setSystem(int system);
	void calculateBidentateMap();
	void setAllRotations(const std::vector<double> &allRotationsVector);
	std::vector<int> calculateRotationTransform(int rotation);

	std::vector<int> applyRotation(const std::vector<int> & permutation, int iRotation);

	std::vector<int> applyPermutationCoordinates(
		const std::vector<int> & permutation,
		const std::vector<std::string> & atoms,
		const std::vector<int> & bidentateAtomsChosen);

	std::vector<int> applyPermutationBidentates(
		const std::vector<int> & permutation,
		const std::vector<int> & bidentateAtomsChosen);

	unsigned int factorial(unsigned int n);

	void printCauchyNotation(std::vector<int> & cauchyList);
	void printCauchyNotation(
		std::string fileName, 
		std::vector<int> & cauchyList);

	void printMoleculeFast(std::vector<CoordXYZ> & mol);

	//no division between rotations - fixed format number nRotations
	size_t nRotations;
	size_t systemSize;
	void writeCauchyRotations(std::string fileName, std::vector< std::vector<int> > & rotPermutations);
	std::vector<int> readCauchyNotations(std::ifstream & openedFile_);

	void molecularFormulaToCauchyCode(
		std::string code,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateAtomsChosen);

	void setBidentateChosen(std::vector<int> & bidentateAtomsChosen);


	AuxMath auxMath_;
};

#endif


