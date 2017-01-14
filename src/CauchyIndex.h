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

private:
	//data rotations
	std::vector<CoordXYZ> mol0;
	std::vector<CoordXYZ> molBidentate;
	std::vector< std::vector<int> > bidentateMap;
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

	std::vector<int> applyPermutation(
		const std::vector<int> & permutation,
		const std::vector<std::string> & atoms,
		const std::vector<int> & bidentateAtomsChosen);

	unsigned int factorial(unsigned int n);

	void printCauchyNotation(std::vector<int> & cauchyList);

	AuxMath auxMath_;
};

#endif