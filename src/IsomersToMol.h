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

	void printAllMol(std::string fileName);

	std::vector<std::string> readAllPermutations(
		std::string fileName, 
		std::vector<int> &atomTypes,
		std::vector<int> &bidentateChosen);

	void printSingleMol(
		std::vector<int> &permutation,
		std::vector<int> &atomTypes,
		std::vector<int> &bidentateChosen,
		std::string name);

private:

	// Functions
	void readAtomTypesAndBidentateChosenFile(
		std::ifstream & file_,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen,
		int systemSize,
		int nBidentates);
	std::string permutationToString0Correction(std::vector<int> &permutation);
	std::string permutationToString(std::vector<int> &permutation);
	std::vector<int> readCauchyNotationsEnantiomers(std::ifstream & openendFile_);
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

	void setParameters(int coordination);

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
