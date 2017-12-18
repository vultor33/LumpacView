#ifndef MARQUESENANTIOMERS_H
#define MARQUESENANTIOMERS_H

#include <vector>
#include <string>

#include "Coordstructs.h"

class MarquesEnantiomers
{
public:
	MarquesEnantiomers();

	~MarquesEnantiomers();

	double marquesRmsd(
		std::vector<CoordXYZ> &mol1,
		std::vector<CoordXYZ> &mol2);

	void calculateMarquesEnantiomers(
		std::string fileXyz1,
		std::string fileXyz2);

	void setSeed(int seed);

private:
	double assignStruct(
		std::vector<CoordXYZ> &mo1,
		std::vector<CoordXYZ> &mol2);

	double setMass(std::string atomLabel);

	void translateToCenterOfMass(std::vector<CoordXYZ> &mol);

	std::vector<CoordXYZ> readXyz(std::string xyzFile);

	void eulerRotation(std::vector<CoordXYZ> &mol);

	double rmsFitConnect(
		std::vector<CoordXYZ> &mol1,
		std::vector<int> &connect1,
		std::vector<CoordXYZ> &mol2,
		std::vector<int> &connect2);

	std::vector< std::vector<double> > distanceMatrix(
		std::vector<CoordXYZ> &mol1,
		std::vector<CoordXYZ> &mol2);

	void jacobiDiagonalization(
		std::vector< std::vector<double> > &entryMatrix,
		std::vector< std::vector<double> > &eigenVectors,
		std::vector<double> &eigenvalues);

	void maxMinBetweenTwoNumber(const int n, const int m, int &iMax, int &iMin);

	void quatFit(
		const std::vector<CoordXYZ> &mol1,
		const std::vector<int> & connect1,
		std::vector<CoordXYZ> &mol2,
		const std::vector<int> & connect2);

	void printCoordXyz(std::vector<CoordXYZ> &mol);

	//paramters
	double tolRmsFit;
	double elementDiff;

	//assign
	double tolIdentificalStructRmsd;
	int nCycles;
	int minRepeat;	
	int kCount;
};


#endif
