#ifndef CONTROLMOPAC_H
#define CONTROLMOPAC_H

#include <vector>
#include <string>

#include "Coordstructs.h"
#include "ReadInput.h"
#include "ReadMopac.h"

class ControlMopac
{
public:
	ControlMopac();
	~ControlMopac();

	bool findLigandsConfiguration(ReadInput &ri_);

	
private:
	std::string projectName;
	std::string metalName;
	std::string mopacHeader;
	std::string mopacFreq;
	std::string metalParams;
	std::vector< Ligand > ligands;
	double frequency;

	void optimization(ReadInput &ri_);
	void forceCalculation(Ligand &optimizedMol_, ReadInput &ri_);
	void buildMopacInput(vector<Ligand> &ligands, string execOption);





};

#endif





