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

	void formatToSymmetryAndFiles(int geoCode);

	void identifyAll();

	void generateAllIsomersMol2Files();

	//testing
	void identifyOne();

	void findAllGroupPoint(int geoCode);

private:

	std::string getResponseName(int size);

	std::string getCountingPath(int geoCode);

	std::string getResultsPath(int geoCode);

	std::string getResultsPathLinux(int geoCode);
};



#endif
