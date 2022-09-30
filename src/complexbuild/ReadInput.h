#ifndef READINPUT_H
#define READINPUT_H

#include <string>
#include <vector>

#include "Coordstructs.h"
#include "Ligand.h"

class ReadInput
{
public:
	ReadInput();
	~ReadInput();

	void readLumpacViewInput();

	void setLumpacViewInputWithoutFile(std::vector< std::string > & inputInformations);

	void rePrintInput();

	void buildLumpacViewFromNames(std::vector< std::string > & fileNames);

	std::vector<Ligand> allLigands;

	inline void setProperties(std::vector<std::string> options_in, std::string mopacExecPath_in)
	{
		options = options_in;
		mopacExecPath = mopacExecPath_in;
	}

	inline std::vector<std::string> getOptions() { return options; }

	inline std::string getMopacExecPath() { return mopacExecPath; }

private:
	std::string inputName;
	//data
	std::vector<std::string> options;
	std::vector< std::string > ligandFileName;
	std::string mopacExecPath;

	std::string buildProjectName(std::string metalName);
	Ligand readConfigurations(std::string inputName);
	void readAllLigands();
	bool failed;

};

#endif

//const std::string mopacHeader = "RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
//const std::string mopacFreq = "RM1 PRECISE NOINTER XYZ T=10D AUX THERMO FORCE + \n NOLOG GEO-OK SCFCRT=1.D-10";
