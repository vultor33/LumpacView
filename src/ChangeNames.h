#ifndef CHANGENAMES_H
#define CHANGENAMES_H

#include <vector>
#include <string>
#include <fstream>

struct vultorGroup
{
	int chiralProb;
	int chiralN;
	int achiralProb;
	int achiralN;
	std::string blockName;
};

class ChangeNames
{
public:
	ChangeNames();
	~ChangeNames();

	void changeNameOfFiles(
		std::string name,
		int geoCode,
		std::string pathRead,
		std::string pathWrite);

	void createNewCounting(
		int geoCode,
		std::string pathRead,
		std::string responseName);

	void generateOrderingGroupPoint(
		std::string fileName,
		std::vector<int> &uniqRcw,
		std::vector<int> &uniqCount,
		std::vector<int> &allNumbers,
		std::vector<std::string> &uniqPgroup);

	void generateOrderingGroupPoint(
		std::string fileName,
		std::vector<int> &uniqRcw,
		std::vector<int> &uniqCount,
		std::vector<int> &allNumbers,
		std::vector<std::string> &uniqPgroup,
		std::string &firstLine,
		std::vector<std::string> & listOfPermutations,
		std::vector<std::string> & listOfChiralities,
		std::vector<int> & listOfRcw,
		std::vector<std::string> & listOfPGroup);

	int symmetryNumberFromPointGroup(std::string pGroup);

private:

	void rceOrderingCriteria(
		std::vector<int> &rcw,
		std::vector<int> &rce,
		std::vector<int> &count,
		std::vector<std::string> &pGroup,
		std::vector<std::string> &setsLetter);

	std::vector<vultorGroup> setVultorGroup(
		std::vector<int> &probs,
		std::vector<int> &number,
		std::vector<bool> &chiral);

	std::string takeLetter(int nGroup);

	vultorGroup findVultorGroup(int prob, std::vector<vultorGroup> &group);

	void calculateVultorGroup(
		std::ifstream & isomerFile_,
		std::vector<vultorGroup> & vgroup);

	void changeFormat(
		std::ifstream &isomerFile_,
		std::ofstream &counting,
		std::ofstream &newFile_,
		std::vector<vultorGroup> &vGroup,
		std::string &geomName,
		std::string &newCombinationName,
		int systemSize
		);

	int calculateSystemSize(std::vector< std::vector<int> > & combinationCode);

	std::string generateNewTypeLine(
		std::string pathRead,
		std::string &combination, 
		int systemSize);

};

#endif

