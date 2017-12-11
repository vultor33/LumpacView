#ifndef ENANTIOMERIDENTIFICATION_H
#define ENANTIOMERIDENTIFICATION_H

#include <vector>

#include "Coordstructs.h"
#include "AuxMath.h"

class EnantiomerIdentification
{
public:
	EnantiomerIdentification();
	~EnantiomerIdentification();

	std::vector<int> finalIsomer(
		std::vector<CoordXYZ> &mol0,
		std::vector<int> &atomTypes,
		std::vector<int> &bidentate,
		std::vector<int> &permut1,
		std::vector<int> &permut2,
		std::vector<CoordXYZ> &crystalGeom,
		std::vector<int> &composition,
		std::vector<int> &bidentateCrystal);

private:
	AuxMath auxMath_;


	std::vector<CoordXYZ> rotateVectorToZ(
		int iAtom,
		std::vector<CoordXYZ> &mol);

	int chooseBestCompareType(std::vector<int> &atomTypes1);

	//first element is of atomTypes1, the next are all corresponding atomTypes2 possible positions
	std::vector<int> mapNeededPermutation(
		int minCount,
		std::vector<int> &atomTypes1,
		std::vector<int> &atomTypes2);

	double compareEnantiomersThroughRotation(
		int chosen1,
		int chosen2,
		std::vector<CoordXYZ> &coord1,
		std::vector<int> &atomTypes1,
		std::vector<int> &bidentate1,
		std::vector<CoordXYZ> &coord2,
		std::vector<int> &atomTypes2,
		std::vector<int> &bidentate2);

	std::vector< std::vector<double> > generateDistanceMatrix(
		std::vector<CoordXYZ> &coord1,
		std::vector<CoordXYZ> &coord2);

	int searchBidentateMatch(
		std::vector<int> &bidentate,
		int j);

	double totalDist(
		int chosen1,
		int chosen2,
		std::vector<CoordXYZ> &coord1,
		std::vector<CoordXYZ> &coord2,
		std::vector< std::vector<double> > &matrDist,
		std::vector<int> atomTypes1,
		std::vector<int> bidentate1,
		std::vector<int> atomTypes2,
		std::vector<int> bidentate2,
		std::vector<int> &connect);

	void rotateCoordOnChosen(
		int chosen,
		std::vector<CoordXYZ> &coord,
		double angle);

	void applyPermutationBidentates(
		const std::vector<int>& permutation,
		std::vector<int>& bidentateAtomsChosen);

};


#endif
