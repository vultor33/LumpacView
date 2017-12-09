#ifndef IDENTIFYISOMERS_H
#define IDENTIFYISOMERS_H

#include <vector>
#include <string>
#include "Coordstructs.h"

/*
EXEMPLO:
ci234.identifyIsomer(
"SAPR-8-M01M01M02M02M02M02M02M02.csv",
"GAKYUM.xyz");

FORMATO DA GEOMETRIA

E necessario conferir o formato no arquivo de permutacoes.
Se o m01 e tipo 1, entao o monodentado que aparecer tem q ser tipo 1.
Os ligantes bidentados vem no final na ordem direta com quebra de linha.
Tipos: 0-Ca; 1-O; 2-H; 3-Mg; 4-Be; 5-I; 6-Cl; 7-Br; 8-Na; 9-He; 10-N; 11-C;
*/

class IdentifyIsomers
{
public:
	IdentifyIsomers();
	~IdentifyIsomers();

	void coordinatesToPermutation(
		std::vector<CoordXYZ> & mol0,
		std::string permutationsFile,
		std::string coordinatesFile);

	double compareGeometryPermutation(
		std::vector<int> &atomTypes,
		std::vector<int> &permutation,
		std::vector<CoordXYZ> &mol0,
		std::vector<int> &bidentatePermutation,
		std::vector<int> &composition,
		std::vector<CoordXYZ> & outGeometry,
		std::vector<int> &bidentateGeometry);

	void compareTwoPermutations(
		std::vector<int> &atomTypes,
		std::string &permutation1,
		std::string &permutation2,
		std::vector<CoordXYZ> &mol0);

private:
	void compareTwoPermutations(
		std::vector<int> &atomTypes,
		std::vector<int> &permutation1,
		std::vector<int> &permutation2,
		std::vector<CoordXYZ> &mol0);

	int indexofSmallestElement(std::vector<double> &vector);

	void insertSortedType(
		int type,
		std::vector<int> & reducedTypes,
		std::vector<double> & distances,
		std::vector<int> & newTypes,
		std::vector<double> & newDistances);

	double compareTwoAtoms(
		std::vector<int> aT1,
		std::vector<double> aD1,
		std::vector<int> aT2,
		std::vector<double> aD2);


	void compareIsomers(
		std::vector<int> & atomTypes1,
		std::vector<int> & bidentates1,
		std::vector< std::vector<int> > & allSortedTypes1,
		std::vector< std::vector<double> > & allSortedDistances1,
		std::vector<int> & atomTypes2,
		std::vector<int> & bidentates2,
		std::vector< std::vector<int> > & allSortedTypes2,
		std::vector< std::vector<double> > & allSortedDistances2,
		std::vector<int> &connections,
		std::vector<double> &difference);


	void sortAllDistances(
		std::vector<int> &atomTypes,
		std::vector<CoordXYZ> &mol0,
		std::vector< std::vector<int> > &allSortedTypes,
		std::vector< std::vector<double> > &allSortedDistances);

	void calculateDistancesK(
		int k,
		std::vector<CoordXYZ> &mol0,
		std::vector<double> &distances);

	void sortDistancesK(
		int k,
		std::vector<int> &atomTypes,
		std::vector<CoordXYZ> &mol0,
		std::vector<int> &sortedTypes,
		std::vector<double> &sortedDistances);

	std::vector<int> stringToPermutation(std::string entryString, size_t size);

	std::vector<CoordXYZ> readGeometry(
		std::string fileName, 
		std::vector<int> &composition,
		int bidSize,
		std::vector<int> &bidentateGeometry);

	void translateToCenterOfMassAndReescale(std::vector<CoordXYZ> &coord);

	int reverseComposition(std::string comp);

	bool compareConnectBidentates(
		std::vector<int> &connect,
		std::vector<int> &bidPermut,
		std::vector<int> &bidGeom);

	int searchBidentateMatch(std::vector<int> &bidentate, int j);

	void applyPermutationBidentates(
		const std::vector<int>& permutation,
		std::vector<int>& bidentateAtomsChosen);

	double differenceConnect(
		std::vector<int> & connect,
		std::vector<int> & atomTypes1,
		std::vector<int> & bidentates1,
		std::vector< std::vector<int> > & allSortedTypes1,
		std::vector< std::vector<double> > & allSortedDistances1,
		std::vector<int> & atomTypes2,
		std::vector<int> & bidentates2,
		std::vector< std::vector<int> > & allSortedTypes2,
		std::vector< std::vector<double> > & allSortedDistances2);


};

#endif
