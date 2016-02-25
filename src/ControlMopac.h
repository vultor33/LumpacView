#ifndef CONTROLMOPAC_H
#define CONTROLMOPAC_H

#include <vector>
#include <string>

#include "Coordstructs.h"
#include "Ligand.h"
#include "ReadMopac.h"

class ControlMopac
{
public:
	ControlMopac();
	~ControlMopac();

	bool optimize(std::vector<Ligand> & ligands);

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





