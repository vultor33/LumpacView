#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include <vector>
#include <string>
#include <fstream>

class Combinations
{
public:
	Combinations();

	~Combinations();

	void doAllCombinations(int nCoordination);

	void findAllIsomersOnCombinations(std::string combinationFile);

private:
	//data
	int coordination;
	std::ofstream printCombinations_;
	std::vector< std::string > elem;
	std::vector<int> coord;
	int printCoord;
	int maxSize;
	std::vector< std::string > allLigandsNames;
	int ligandsSeparationSize;

	//functions
	void clearEqualCombinations(std::string combFile);

	std::vector< std::vector<int> > stringToNumber(std::string entryString);

	std::string codeToString(std::vector< std::vector<int> > & codeLine);

	bool compareToAll(std::vector< std::vector< std::vector<int> > > & allCodes, std::vector< std::vector<int> > &actualCodes);

	int codeToType(std::string code);

	void printInputWithCode(std::string code);

	void codeToLigands(
		std::string code,
		std::vector< std::string > & ligandNames,
		std::vector<int> & denticity);

	void addEqual(
		int codeNumber,
		std::vector<int> & typeCode1,
		std::vector<int> & typeCode2,
		std::vector<int> & typeCode3
		);

	void addDifferent(
		int codeNumber,
		std::vector<int> & typeCode1,
		std::vector<int> & typeCode2,
		std::vector<int> & typeCode3
		);

	void sumUpNameAndPrint(
		std::vector< std::string > & name,
		std::vector< int > & totalCoord,
		int position,
		std::string newLigand,
		int coordination);

	int angleBidentateCut(size_t coordination);

};


#endif
