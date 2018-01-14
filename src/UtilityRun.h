#ifndef UTILITYRUN_H
#define UTILITYRUN_H

#include <string>

class UtilityRun
{
public:
	UtilityRun();

	~UtilityRun();

	void renameAtomTypes(std::string responseName);

	void formatIsomersFiles();

	void identifyAll();

	void generateAllIsomersMol2Files();

	//testing
	void identifyOne();


private:

	std::string getCountingPath(int geoCode);

	std::string getResultsPath(int geoCode);
};



#endif