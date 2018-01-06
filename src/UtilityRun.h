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

	void identifyOne();

	void identifyAll();

private:

	std::string getCountingPath(int geoCode);


};



#endif