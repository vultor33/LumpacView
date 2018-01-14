#ifndef ISOMERSTOMOL_H
#define ISOMERSTOMOL_H

#include <vector>
#include <string>
#include <fstream>

#include "Coordstructs.h"

class IsomersToMol
{
public:
	IsomersToMol();

	~IsomersToMol();

	void printAllMolFromSpecifiedGeometry(
		int geoCode,
		std::string pathRead,
		std::string responseName);

	std::vector<std::string> readAllPermutations(
		std::string fileName, 
		std::string fileFolder,
		std::vector<int> &atomTypes,
		std::vector<int> &bidentateChosen);

	std::vector<int> findEnantiomerPair(
		std::string fileName,
		std::string fileFolder,
		std::vector<int> guessPermutation,
		std::vector<std::string> &pairCodes);

	void printSingleMol(
		std::vector<int> &permutation,
		std::vector<int> &atomTypes,
		std::vector<int> &bidentateChosen,
		std::string name);

	void printCoordXyz(std::vector<CoordXYZ> &coord);

	std::string getAtomLabelI(int I);
	
	/* fredapagar
	void readAtomTypesAndBidentateChosenFile(
		std::ifstream & file_,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen,
		int systemSize,
		int nBidentates);
	*/

	std::vector<int> readCauchyNotationsEnantiomers(std::ifstream & openendFile_);

private:
	void printAllMol(
		std::string fileName,
		std::string filePath,
		int geoCode);

	// Functions
	std::string permutationToString0Correction(std::vector<int> &permutation);
	std::string permutationToString(std::vector<int> &permutation);
	std::vector<int> readCauchyNotationsEnantiomersAndTakeCode(
		std::ifstream & openendFile_,
		std::vector<std::string> &permutCodes);

	void printMoleculeMolFormat(
		std::vector<int> & permutation,
		std::vector<int> & atomTypes,
		const std::vector<int> & bidentateAtomsChosen,
		std::ofstream & printFile_);
	std::vector<int> applyPermutationCoordinates(
		const std::vector<int> & permutation,
		const std::vector<std::string> & atoms,
		const std::vector<int> & bidentateAtomsChosen);

	// read composition
	/* fredapagar
	int stringToNumber(std::string entryString, int &nBidentates);
	int codeToType(std::string code);
	void addEqual(
		int codeNumber,
		std::vector<int> & typeCode1,
		std::vector<int> & typeCode2,
		std::vector<int> & typeCode3);
	void addDifferent(
		int codeNumber,
		std::vector<int> & typeCode1,
		std::vector<int> & typeCode2,
		std::vector<int> & typeCode3);
	*/

	// raswin parameters
	void setParameters(int coordination);

	// molekel parameters
	void setParameters(int coordination, int geoCode);

	std::vector<std::string> atomLabels;

	std::vector<CoordXYZ> complex;

};

#endif


/*
cout << "Input name:  ";
getline(cin, fileName);

ifstream fileIsomers_((fileName).c_str());
if (!fileIsomers_.good())
{
cout << "Couldn't find:   " << fileName << endl
<< "Press enter to exit" << endl;
string dummy;
getline(cin, dummy);
exit(1);
}
*/
