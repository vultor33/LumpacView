#ifndef CONTROLMOPAC_H
#define CONTROLMOPAC_H

#include <vector>
#include <string>

#include "Coordstructs.h"
#include "ReadMopac.h"

class ControlMopac
{
public:
	ControlMopac(
		std::string projectName_in,
		std::string metalName_in,
		std::string mopacHeader_in,
		std::string mopacFreq_in,
		std::string metalParams_in
		);
	~ControlMopac();

	bool optimize(std::vector<CoordXYZ> & allAtoms);

private:
	std::string projectName;
	std::string metalName;
	std::string mopacHeader;
	std::string mopacFreq;
	std::string metalParams;
	std::string mopacExecPath = "M2009_Ln_Orbitals.exe  ";
	std::vector< Ligand > ligands;
	double frequency;

	double forceCalculation(std::vector<CoordXYZ> &optimizedAtoms);
	void buildMopacInput(std::vector<CoordXYZ> &allAtoms, std::string execOption);
	
};

#endif

