#include "ReadMopac.h"

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "Coordstructs.h"
#include "Ligand.h"

using namespace std;

ReadMopac::ReadMopac(){}

ReadMopac::~ReadMopac(){}

bool ReadMopac::readOutput(string outName, Ligand &mol_)
{
	string name = outName + ".out";
	ifstream in_;
	in_.open(name.c_str());

	vector<CoordXYZ> newCoord;
	bool first = false;
	bool second = false;
	bool readCoord = false;
	int natoms;
	string auxline, auxCartesian, auxCoordinates, dum1, dum2;
	while (getline(in_, auxline))
	{
		if (!second)
		{
			stringstream line;
			line << auxline;
			line >> auxCartesian >> auxCoordinates;
			if ((auxCartesian == "CARTESIAN") &&
				(auxCoordinates == "COORDINATES")
				)
			{
				if (first)
				{
					second = true;
					getline(in_, auxline);
					getline(in_, auxline);
				}
				else
				{
					auxCartesian = "";
					auxCoordinates = "";
					first = true;
				}
			}
		}
		else
		{
			if (readCoord)
			{
				if (auxline == "")
				{
					break;
				}
				stringstream line;
				line << auxline;
				CoordXYZ auxcoord;
				line >> natoms >> auxcoord.atomlabel>> auxcoord.x
					>> auxcoord.y >> auxcoord.z;

				newCoord.push_back(auxcoord);
				//std::cout << "nat " << natoms << " lab " << auxcoord.atomlabel
				//<< "  " << auxcoord.x << "   " << auxcoord.y
				//<< "  " << auxcoord.z << endl;
			}
			else
			{
				readCoord = true;
			}
		}
	}
	in_.close();

	mol_.setNewCoordinates(newCoord);

	return (newCoord.size() != 0);
}



double ReadMopac::readFrequency(string freqName)
{
	string name = freqName + ".out";
	ifstream in_;
	in_.open(name.c_str());

	string auxline, dummy, dummy2, dummy3;

	while (getline(in_, auxline))
	{
		if (auxline == "")
			continue;

		stringstream line0;
		line0 << auxline;
		line0 >> dummy >> dummy2 >> dummy3;

		if ((dummy == vibrationFlag1)&&
			(dummy2 == vibrationFlag2)&&
			(dummy3 == vibrationFlag3))
		{
			getline(in_, auxline);
			getline(in_, auxline);
			getline(in_, auxline);
			getline(in_, auxline);
			stringstream line1;
			double auxFrequency;
			line1 << auxline;
			line1 >> dummy >> auxFrequency;

			in_.close();
			return auxFrequency;
		}
	}
	in_.close();

	return -100;
}



double ReadMopac::readTotalEnergy(string outName)
{
	string name = outName + ".out";
	ifstream in_;
	in_.open(name.c_str());

	string auxline, dummy, dummy2, dummy3;
	double energy;

	while (getline(in_, auxline))
	{
		if (auxline == "")
			continue;

		stringstream line0;
		line0 << auxline;
		line0 >> dummy >> dummy2 >> dummy3;

		if ((dummy == energyFlag1) &&
			(dummy2 == energyFlag2) &&
			(dummy3 == energyFlag3))
		{
			line0 >> energy;
			return energy;
		}
	}
	in_.close();

	return 1.0e0;
}

