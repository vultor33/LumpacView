#ifndef CAUCHYINDEX_H
#define CAUCHYINDEX_H

#include <vector>
#include <string>

#include "Coordstructs.h"
#include "AuxMath.h"

struct cauchyRotation
{
	std::vector< std::vector<double> > mRot;
};

class CauchyIndex
{
public:
	CauchyIndex();
	
	~CauchyIndex();

	void getAllRotationMatrix();

	std::vector<int> getCauchy(std::string fName);

	void calculateAllIndexes();

	void printMolecule(
		std::vector<int> & permutation,
		std::vector<std::string> & atoms,
		std::vector<int> bidentateAtomsChosen);

private:
	std::vector<CoordXYZ> mol0;

	std::vector<CoordXYZ> molBidentate;

	std::vector<cauchyRotation> allRotations;

	std::vector< std::vector<int> > bidentateMap;

	void calculateBidentateMap();

	void setSystem(int system);

	void setAllRotations(const std::vector<double> &allRotationsVector);

	void printCauchyNotation(std::vector<int> & cauchyList);

	std::vector<int> getCauchy(int rotation);

	std::vector<std::string> bidentateLabels;

	AuxMath auxMath_;


};

#endif