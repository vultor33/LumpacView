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

	void changeNameOfFiles(std::string name);

private:
	std::vector<vultorGroup> setVultorGroup(
		std::vector<int> &probs,
		std::vector<int> &number,
		std::vector<bool> &chiral);

	std::string takeLetter(int nGroup);

	vultorGroup findVultorGroup(int prob, std::vector<vultorGroup> &group);

	std::string sizeToGeometryCode(int size);

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

	std::string generateNewTypeLine(std::string &combination, int systemSize);

};

#endif

