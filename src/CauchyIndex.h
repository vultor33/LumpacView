#ifndef CAUCHYINDEX_H
#define CAUCHYINDEX_H

#include <vector>
#include <string>
#include <fstream>
//#include <unistd.h>
#include <sys/stat.h>

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
		std::vector<int> & atomTypes,
		const std::vector<int> & bidentateAtomsChosen,
		std::ofstream & printFile_);
	void printMolecule(
		std::vector<int> & permutation,
		const std::vector<std::string> & atoms,
		const std::vector<int> & bidentateAtomsChosen,
		std::ofstream & printFile_);

	void printMoleculeFromFile(std::string fileName);

	void printAllMoleculesFromFile(std::string composition);

	void rotationTest(
		std::vector<std::string> &atoms,
		std::vector<int> & bidentateAtomsChosen);

	void generateAllIndependentIsomers();

	void generateAllIndependentIsomersIO();

	void generateAllIndependentIsomersRuntimeRotations();

	void generateAllIndependentIsomersRuntimeRotationsAndReadBlock(std::string blockFileName);

	void generateAllIndependentIsomersWithFlag(
		std::string blockFileName, 
		std::string flagsFile,
		std::string code);

	std::vector<int> zeroPermutation(std::string flagsFile);

	void correctIndependentIsomers();

	void generateAllIndependentIsomers12(std::string blockFileName);

	//0-true ; 1-false ; 2-bidentate problem
	int compareTwoIsomersWithLabelsRotations(
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateAtomsChosen,
		std::vector<int> & permutationIsomer1,
		std::vector<int> & permutationIsomer2
		);

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

	void generateBlockFiles(int n, int kInit, int kFinal);

	void generateRAMBlock(int n, int kInit, int kFinal, std::vector< std::vector<int> > & ramBlock);

	// machineType:  slurm | pc
	void generateSlurmFilesToDeletion(int nSystem, int nProc, std::string machineType);

	void generateSlurmFilesToDeletionFlags(
		        int deletionSystem,
        		int total,
        		int bigBlockSize,
        		int smallBlockSize,
        		std::string compositionFile,
        		std::string rawIsomersFile,
				std::string workingDir,
        		std::string machineType);

	void runall(int blockInit, int blockFinal, std::string machineType, std::string workingDirectory);

	void cleanBlocksAndGenerateIsomers(
		int nProc, 
		int systemSize, 
		std::string composition,
		std::string workingDirectory, 
		std::string machineType);

	void doBlockDeletion(
		int kInit,
		int kFinal);

	void doBlockRAMDeletion(
		int kInit,
		int kFinal,
		int ramInit,
		int ramFinal);


	void doBlockRAMDeletion12(
		int kInit,
        	int kFinal,
        	int ramInit,
       		int ramFinal);

	void generateAtomTypesAndBidentateChosenFile(std::string complexCode);

	void readAtomTypesAndBidentateChosenFile(
		std::string fileName,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen);

	void doBlockDeletionFlags(
		std::string skeletonFile,
		std::string flagsFile,
		int kInit,
		int kFinal,
		int ramInit,
		int ramFinal);

	void enantiomersOrdering();

	void enantiomersOrderingBlock(
		int iPermut,
		int kFinal,
		std::string fileName);

// uma parada que leia os isomeros esqueletos
// ele precisa renumerar. posso soltar so os 120 mesmo.
// os flags precisam ser os mesmos sempre - vou ter que gerar um arquivo na pasta para todos lerem.



private:
	//data rotations
	std::vector<CoordXYZ> mol0;
	std::vector<CoordXYZ> molBidentate; // all bidentate possible positions marked with atoms.
	std::vector< std::vector<int> > bidentateMap; // between i and j we have k. K position is the same shown on molBidentate.
	std::vector<cauchyRotation> allRotationsMatrix;
	std::vector< std::vector<int> > allRotationTransforms;
	std::vector<std::string> bidentateLabels;
	std::vector<std::string> atomLabels;
	std::vector<int> reflectionOperation;
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

	std::vector<int> applyZAxisReflection(std::vector<int> permutation);

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
	std::vector< std::vector<int> > readCauchyNotationsRAMBlock(std::ifstream & openedFile_, int kInitial, int kFinal);

	std::string permutationToString(std::vector<int> permutation); 

	void molecularFormulaToCauchyCode(
		std::string code,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateAtomsChosen);

	void setBidentateChosen(std::vector<int> & bidentateAtomsChosen);

	inline bool exist_file (const std::string& name)
	{
		struct stat buffer;
  		return (stat (name.c_str(), &buffer) == 0);
	}

	AuxMath auxMath_;
};

#endif


