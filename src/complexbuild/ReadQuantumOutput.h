#ifndef READQUANTUMOUTPUT_H
#define READQUANTUMOUTPUT_H

#include <string>
#include <vector>
#include <fstream>

#include "Coordstructs.h"

class ReadQuantumOutput
{
public:
	ReadQuantumOutput(std::string type_in);
	~ReadQuantumOutput();

	 void readOutput(std::string outName);

	 //get values
	 std::vector<CoordXYZ> getCoordinates();

	 inline double getEnergy() { return energy; }

	 inline double getIonizationPotential() { return ionizationPotential; }

	 std::vector<double> getDipole();

	 inline double getFirstFrequency(){ return firstFrequency; }

	 //set options
	 void activateDeactivateReadings(std::string activateOption, bool activate);

private:
	std::string type;

	//GENERIC DATA
	std::vector<CoordXYZ> coordinates;
	std::vector<double> dipole;
	double energy;
	double ionizationPotential;
	double firstFrequency;
	bool coordinatesActivation;
	bool energyActivation;
	bool ionizationActivation;
	bool dipoleActivation;
	bool frequencyActivation;

	//END OF GENERERIC DATA

	//MOPAC FLAGS
	std::string vibrationFlag1;
	std::string vibrationFlag2;
	std::string vibrationFlag3;
	std::string energyFlag1;
	std::string energyFlag2;
	std::string energyFlag3;
	std::string ionizationFlag1;
	std::string ionizationFlag2;
	std::string ionizationFlag3;
	std::string dipoleFlag1;
	std::string dipoleFlag2;
	std::string dipoleFlag3;
	std::string dipoleFlag4;
	std::string dipoleFlag5;
	std::string coordFlag1;
	std::string coordFlag2;
	bool secondCartesian;
	//END OF MOPAC FLAGS

	//GAMESS FLAGS
	std::string gamessXyzCoordinatesFlag;
	std::string gamessEnergyFlag;
	std::string gamessGradientFlag;
	std::string gamessIonizationFlag;
	std::string gamessDipoleFlag;
	std::string gamessFrequency;
	int gamessIonizationPos;
	bool stopReadingFrequency;
	double readnDoubles(std::string auxline, int nEigens);
	//END OF GAMESS FLAGS

	// Control action - how to read.
	std::vector<CoordXYZ> readCoordinates(std::ifstream & quantumOut_);
	void readEnergy(std::string auxline);
	void readIonization(std::string auxline, std::ifstream & quantumOut_);
	void readDipole(std::ifstream & quantumOut_);
	void readFrequency(std::ifstream & quantumOut_);

	// Control flags - when to read.
	bool haveToReadCoordinates(std::string auxline);
	bool haveToReadEnergy(std::string auxline);
	bool haveToReadIonization(std::string auxline);
	bool haveToReadDipole(std::string auxline);
	bool haveToReadFrequency(std::string auxline);

};

#endif

