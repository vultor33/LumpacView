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

	// {SAPR-8 [Ma2b2c2] [1 2 3 4 5 6 7 8] Aa} --> format
	vector<int> readCauchyNotationsEnantiomers(
		ifstream & openendFile_,
		int size);

	std::string codeToString(std::vector < std::vector<int> > & codeLine);

	std::string newCodeToString(std::vector < std::vector<int> > & codeLine);

	// (m, B e C)  versao do allMolecular
	std::vector< std::vector<int> > stringToNumber(std::string entryString);

	// {a (AA) (AB)}  versao do isomers to mol
	int CompositionToNumbers(
		std::string entryString,
		int & nBidentates);


private:
	std::vector< std::string > elem;
	std::vector< std::string > elemNew;
	
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



};


#endif