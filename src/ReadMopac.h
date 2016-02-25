#ifndef READMOPAC_H
#define READMOPAC_H

#include <string>

#include "Ligand.h"

class ReadMopac
{
public:
	ReadMopac();
	~ReadMopac();

	bool readOutput(std::string outName, Ligand &mol_);

	double readFrequency(std::string freqName);

	double readTotalEnergy(std::string outName);

//	inline double getFirstFrequency(){ return firstFrequency; }

private:
	const std::string vibrationFlag1 = "DESCRIPTION";
	const std::string vibrationFlag2 = "OF";
	const std::string vibrationFlag3 = "VIBRATIONS";
	const std::string energyFlag1 = "TOTAL";
	const std::string energyFlag2 = "ENERGY";
	const std::string energyFlag3 = "=";
};

#endif

