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

	void rePrintInput();

private:
	const std::string inputName = "LumpacViewInput.txt";
	//data
	std::string metalName;
	std::string metalParams;
	std::vector< std::string > ligandFileName;
	std::string projectName;
	std::vector<Ligand> allLigands;

	void readLumpacViewInput();
	void buildProjectName();
	Ligand readConfigurations(std::string inputName);
	void readAllLigands();
	bool failed;

};

#endif

//const std::string mopacHeader = "RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
//const std::string mopacFreq = "RM1 PRECISE NOINTER XYZ T=10D AUX THERMO FORCE + \n NOLOG GEO-OK SCFCRT=1.D-10";
