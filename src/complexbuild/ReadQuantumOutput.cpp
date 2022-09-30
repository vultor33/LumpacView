#include "ReadQuantumOutput.h"

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdlib.h>

#include "Coordstructs.h"

using namespace std;

ReadQuantumOutput::ReadQuantumOutput(string type_in)
{
	type = type_in;
	coordinatesActivation = true;
	energyActivation = true;
	ionizationActivation = false;
	dipoleActivation = false;
	frequencyActivation = true;
	if ((type == "mopac") || (type == "mopac2009"))
	{
		energyFlag1 = "TOTAL";
		energyFlag2 = "ENERGY";
		energyFlag3 = "=";
		dipoleFlag1 = "DIPOLE";
		dipoleFlag2 = "X";
		dipoleFlag3 = "Y";
		dipoleFlag4 = "Z";
		dipoleFlag5 = "TOTAL";
		ionizationFlag1 = "IONIZATION";
		ionizationFlag2 = "POTENTIAL";
		ionizationFlag3 = "=";
		coordFlag1 = "CARTESIAN";
		coordFlag2 = "COORDINATES";
		vibrationFlag1 = "DESCRIPTION";
		vibrationFlag2 = "OF";
		vibrationFlag3 = "VIBRATIONS";
		secondCartesian = false;
	}
	else if (type == "gamess")
	{
		gamessXyzCoordinatesFlag = "COORDINATES OF ALL ATOMS ARE (ANGS)";
		gamessEnergyFlag = "TOTAL ENERGY =";
		gamessGradientFlag = "UNITS ARE HARTREE/BOHR";
		gamessIonizationFlag = "EIGENVECTORS";
		gamessDipoleFlag = "ELECTROSTATIC MOMENTS";
		gamessFrequency = "FREQUENCY:";
		stopReadingFrequency = false;
	}
}

ReadQuantumOutput::~ReadQuantumOutput(){}

vector<CoordXYZ> ReadQuantumOutput::getCoordinates()
{
	return coordinates;
}

vector<double> ReadQuantumOutput::getDipole()
{
	return dipole;
}

void ReadQuantumOutput::activateDeactivateReadings(std::string activateOption, bool activate)
{
	if (activateOption == "coordinates")
		coordinatesActivation = activate;
	else if (activateOption == "energy")
		energyActivation = activate;
	else if (activateOption == "ionization")
		ionizationActivation = activate;
	else if (activateOption == "dipole")
		dipoleActivation = activate;
	else if (activateOption == "frequency")
		frequencyActivation = activate;
	else
	{
		cout << "ERROR ON: ReadQuantumOutput::activateDeactivateReadings" << endl
			<<  "--> option not found" << endl;
		exit(1);
	}
}

void ReadQuantumOutput::readOutput(string outName)
{
	firstFrequency = -1.0e99;
	energy = 0.0e0;
	ionizationPotential = 0.0e0;
	string name = outName + ".out";
	ifstream quantumOut_;
	quantumOut_.open(name.c_str());

	string auxline;
	while (getline(quantumOut_, auxline))
	{
		if (auxline == "")
			continue;

		if (haveToReadCoordinates(auxline) && coordinatesActivation)
		{
			coordinates = readCoordinates(quantumOut_);
		}
		else if (haveToReadEnergy(auxline) && energyActivation)
		{
			readEnergy(auxline);
		}
		else if (haveToReadIonization(auxline) && ionizationActivation)
		{
			readIonization(auxline, quantumOut_);
		}
		else if (haveToReadDipole(auxline) && dipoleActivation)
		{
			readDipole(quantumOut_);
		}
		else if (haveToReadFrequency(auxline) && frequencyActivation)
		{
			readFrequency(quantumOut_);
		}
	}

	quantumOut_.close();
}

vector<CoordXYZ> ReadQuantumOutput::readCoordinates(ifstream & quantumOut_)
{
	string auxline;
	vector<CoordXYZ> newCoord;
	if ((type == "mopac") || (type == "mopac2009"))
	{
		int dummy;
		getline(quantumOut_, auxline);
		if (type == "mopac2009")
		{
			getline(quantumOut_, auxline);
			getline(quantumOut_, auxline);
		}

		while (getline(quantumOut_, auxline))
		{
			if (auxline == "")
				break;

			stringstream line;
			line << auxline;
			CoordXYZ auxcoord;
			line >> dummy >> auxcoord.atomlabel >> auxcoord.x
				>> auxcoord.y >> auxcoord.z;
			newCoord.push_back(auxcoord);
		}
	}
	else if (type == "gamess")
	{
		string auxline;
		getline(quantumOut_, auxline);//dummy line
		getline(quantumOut_, auxline);//dummy line
		string dummy;
		while (getline(quantumOut_, auxline))
		{
			if (auxline == "")
				break;

			stringstream line;
			line << auxline;
			CoordXYZ auxcoord;
			line >> auxcoord.atomlabel >> dummy >> auxcoord.x
				>> auxcoord.y >> auxcoord.z;
			newCoord.push_back(auxcoord);
		}
	}

	return newCoord;
}


void ReadQuantumOutput::readEnergy(string auxline)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		stringstream line;
		string aux1, aux2, aux3;
		line << auxline;
		line >> aux1 >> aux2 >> aux3 >> energy;
	}
	if (type == "gamess")
	{
		stringstream line;
		line << auxline;
		string aux1;
		string aux2;
		string aux3;
		line >> aux1 >> aux2 >> aux3 >> energy;
	}
}

void ReadQuantumOutput::readIonization(string auxline, ifstream & quantumOut_)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		stringstream line;
		string aux1, aux2, aux3;
		line << auxline;
		line >> aux1 >> aux2 >> aux3 >> ionizationPotential;
	}
	if (type == "gamess")
	{
		string auxline;
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);

		bool auxRead = true;
		int k = gamessIonizationPos;
		while (getline(quantumOut_, auxline))
		{
			if (auxRead)
			{
				if (k < 6)
				{
					ionizationPotential = readnDoubles(auxline, k);
				}
				else
				{
					k -= 5;
					auxRead = false;
				}
			}
			else
			{
				if (auxline.size() == 0)
				{
					getline(quantumOut_, auxline);
					auxRead = true;
				}
			}
		}

	}
}

void ReadQuantumOutput::readDipole(ifstream & quantumOut_)
{
	string auxline;
	if ((type == "mopac") || (type == "mopac2009"))
	{
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		stringstream line;
		string aux1;
		line << auxline;
		double x, y, z, total;
		line >> aux1 >> x >> y >> z >> total;
		dipole.resize(4);
		dipole[0] = x;
		dipole[1] = y;
		dipole[2] = z;
		dipole[3] = total;
	}
	if (type == "gamess")
	{
		dipole.resize(4);
		double dx, dy, dz, d;
		string auxline;
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		stringstream convert;
		convert << auxline;
		convert >> dx >> dy >> dz >> d;
		dipole[0] = dx;
		dipole[1] = dy;
		dipole[2] = dz;
		dipole[3] = d;
	}
}

void ReadQuantumOutput::readFrequency(ifstream & quantumOut_)
{
	string auxline;
	if ((type == "mopac") || (type == "mopac2009"))
	{
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		getline(quantumOut_, auxline);
		stringstream line;
		string aux1;
		line << auxline;
		line >> aux1 >> firstFrequency;
	}
	if (type == "gamess")
	{
		while(getline(quantumOut_, auxline))
		{
			if (auxline.find(gamessFrequency) != string::npos)
			{
				string dum1, dum2;
				stringstream convert;
				convert << auxline;
				convert >> dum1 >> dum2 >> firstFrequency;
				break;
			}
		}
	}
}


bool ReadQuantumOutput::haveToReadCoordinates(string auxline)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		string aux1, aux2;
		stringstream line;
		line << auxline;
		line >> aux1 >> aux2;
		bool found = (aux1 == coordFlag1) && (aux2 == coordFlag2);
		if (found)
		{
			if (!secondCartesian)
			{
				secondCartesian = true;
				return false;
			}
			else
				return true;
		}
		else
			return false;
	}
	else if (type == "gamess")
	{
		return auxline.find(gamessXyzCoordinatesFlag) != string::npos;
	}
	
	return false;
}

bool ReadQuantumOutput::haveToReadEnergy(string auxline)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		string aux1, aux2, aux3;
		stringstream line;
		line << auxline;
		line >> aux1 >> aux2 >> aux3;
		return (aux1 == energyFlag1)
			&& (aux2 == energyFlag2)
			&& (aux3 == energyFlag3);
	}
	else if (type == "gamess")
	{
		return auxline.find(gamessEnergyFlag) != string::npos;
	}

	return false;
}

bool ReadQuantumOutput::haveToReadIonization(string auxline)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		string aux1, aux2, aux3;
		stringstream line;
		line << auxline;
		line >> aux1 >> aux2 >> aux3;
		return (aux1 == ionizationFlag1)
			&& (aux2 == ionizationFlag2)
			&& (aux3 == ionizationFlag3);
	}
	else if (type == "gamess")
	{
		return auxline.find(gamessIonizationFlag) != string::npos;
	}

	return false;
}

bool ReadQuantumOutput::haveToReadDipole(string auxline)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		string aux1, aux2, aux3, aux4, aux5;
		stringstream line;
		line << auxline;
		line >> aux1 >> aux2 >> aux3 >> aux4 >> aux5;
		return (aux1 == dipoleFlag1)
			&& (aux2 == dipoleFlag2)
			&& (aux3 == dipoleFlag3)
			&& (aux4 == dipoleFlag4)
			&& (aux5 == dipoleFlag5);
	}
	else if (type == "gamess")
	{
		return auxline.find(gamessDipoleFlag) != string::npos;
	}

	return false;
}

bool ReadQuantumOutput::haveToReadFrequency(string auxline)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		string aux1, aux2, aux3;
		stringstream line;
		line << auxline;
		line >> aux1 >> aux2 >> aux3;
		return (aux1 == vibrationFlag1)
			&& (aux2 == vibrationFlag2)
			&& (aux3 == vibrationFlag3);
	}
	else if (type == "gamess")
	{
		if (stopReadingFrequency)
			return false;

		if (auxline.find(gamessFrequency) != string::npos)
		{
			stopReadingFrequency = true;
			return true;
		}	
	}

	return false;
}

double ReadQuantumOutput::readnDoubles(string auxline, int nEigens)
{
	stringstream convert;
	convert << auxline;
	double a1, a2, a3, a4, a5;
	switch (nEigens)
	{
	case 1:
		convert >> a1;
		return a1;
		break;

	case 2:
		convert >> a1 >> a2;
		return a2;
		break;

	case 3:
		convert >> a1 >> a2 >> a3;
		return a3;
		break;

	case 4:
		convert >> a1 >> a2 >> a3 >> a4;
		return a4;
		break;

	case 5:
		convert >> a1 >> a2 >> a3 >> a4 >> a5;
		return a5;
		break;

	default:
		cout << "Error in readnDoubles - contact developers" << endl;
		exit(1);
	}
}



// counting time
/*
#include <time.h>
time_t rawtime;
struct tm * timeinfo;
time(&rawtime);
timeinfo = localtime(&rawtime);
printf("Current local time and date: %s", asctime(timeinfo));
*/
