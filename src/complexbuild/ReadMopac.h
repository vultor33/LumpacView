#ifndef READMOPAC_H
#define READMOPAC_H

#include <string>

#include "Ligand.h"

class ReadMopac
{
public:
	ReadMopac();
	~ReadMopac();

	std::vector<CoordXYZ> readOutput(std::string outName);

	double readFrequency(std::string freqName);

	double readTotalEnergy(std::string outName);

//	inline double getFirstFrequency(){ return firstFrequency; }

private:
	std::string vibrationFlag1;
	std::string vibrationFlag2;
	std::string vibrationFlag3;
	std::string energyFlag1;
	std::string energyFlag2;
	std::string energyFlag3;
};

#endif

