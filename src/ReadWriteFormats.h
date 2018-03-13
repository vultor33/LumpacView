#ifndef READWRITEFORMATS_H
#define READWRITEFORMATS_H

#include <fstream>
#include <vector>
#include <string>

class ReadWriteFormats
{
public:
	ReadWriteFormats();

	~ReadWriteFormats();

	// number / bidentares   --> format
	void readAtomTypesAndBidentateChosenFile(
		std::ifstream & file_,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen,
		int systemSize,
		int nBidentates);

	// labels / bidentares   --> format
	void readAtomTypesAndBidentateChosenFileWithLabels(
		std::ifstream & file_,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen,
		int systemSize,
		int nBidentates);

	/*
	// letters / bidentares   --> format
	void readAtomTypesAndBidentateChosenFileWithLetters(
		std::ifstream & file_,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen,
		int systemSize,
		int nBidentates);
		*/

	// {SAPR-8 [Ma2b2c2] [1 2 3 4 5 6 7 8] Aa} --> format
	std::vector<int> readCauchyNotationsEnantiomers(
		std::ifstream & openendFile_,
		int size);

	std::string codeToString(std::vector < std::vector<int> > & codeLine);

	std::string newCodeToString(std::vector < std::vector<int> > & codeLine);

	// (m, B e C)  versao do allMolecular
	std::vector< std::vector<int> > compositionToNumberOld(std::string entryString);

	// {a, (AA), (AB)}  versao do isomers to mol
	int compositionToNumbers(
		std::string entryString,
		int & nBidentates);

	// numbers to letters
	std::string typeLineToLetters(
		std::string typeLine,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen);

	// Reverse numbers to letters
	void typeLineToNumberCodes(
		std::string typeLine,
		std::vector<int> & atomTypes,
		std::vector<int> & bidentateChosen);

	std::vector<int> readCauchyNotation(
		std::ifstream & openendFile_,
		int size,
		std::string & vultorGroup,
		int & rcw,
		std::string & allCode);

	void takeRcwVgroupPointGroup(
		std::string line,
		int & rcw,
		std::string & vGroup,
		std::string & pGroup);

	void takeRcwVgroupPointGroupNca(
		std::string line,
		int & rcw,
		int & chiral,
		int & achiral,
		std::string & vGroup,
		std::string & pGroup);

	void symmetryGroupOrdering(
		std::vector<int> &uniqRcw,
		std::vector<std::string> &uniqPgroup,
		std::vector<int> &uniqCount);

	bool hyerarchyOrdering(
		int rcw1,
		int rcw2,
		std::string gPoint1,
		std::string gPoint2);

	bool hyerarchyOrdering(
		std::string gPoint1,
		std::string gPoint2);

	bool hyerarchyOrdering2(
		std::string gPoint1,
		std::string gPoint2);

	std::string includeGroupPoint(
		std::string vCode,
		std::string gPoint);

	std::string takeGroupPoint(std::string vCode);
	
	std::string atomTypesToAtomStrings(
		std::vector<int> atomTypes);

	std::vector<int> atomStringToAtomTypes(
		std::vector<std::string> &atomsStrings);

	void ReplaceAll(std::string & data, char toSearch, std::string addStr);

	// old format
	void takeAllElementsFromCode(
		std::string line,
		int & rcw,
		int & chiral,
		int & achiral,
		std::string & vGroup,
		std::string & pGroup,
		std::string & permut);

	// format {[Ma2(AA)(AB)] OC-6 Cs a B [1 2 4 6 5 3]}
	void takeAllElementsFromCodeNew(
		std::string line,
		int coordination,
		int & rcw,
		std::string & chirality,
		std::string & vGroup,
		std::string & pGroup,
		std::string & setGroup,
		std::vector<int> & permut);



	bool isachiral(std::string pGroup);

private:
	std::vector< std::string > elem;
	std::vector< std::string > elemNew;
	std::vector< std::string > atomLabels;
	
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


	int findCharOnString(
		std::string word,
		char refChar);


};


#endif
