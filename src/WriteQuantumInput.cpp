#include "WriteQuantumInput.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cctype>
#include <stdlib.h>

#include "Coordstructs.h"

using namespace std;

WriteQuantumInput::WriteQuantumInput(
	vector<string> options)
{
	type = options[0];
	projectName = options[1];
	setInputProperties(options);
}

WriteQuantumInput::~WriteQuantumInput() {}

void WriteQuantumInput::setInputProperties(vector<string> &options)
{
	if ((type == "mopac") || (type == "mopac2009"))
	{
		mopacHeader = options[2];
		externalParams = options[3];
		centralMetal = options[4];
	}
	else if (type == "gamess")
	{
		int iEndHeader;
		for (size_t i = 2; i < options.size(); i++)
		{
			if (options[i] == "EndOfHeader")
			{
				iEndHeader = i;
				break;
			}

			gamessHeader.push_back(options[i]);
		}

		int iEndBasis;
		for (size_t i = iEndHeader + 1; i < options.size(); i++)
		{
			if (options[i] == "EndOfBasis")
			{
				iEndBasis = i;
				break;
			}

			gamessAtomBasisFiles.push_back(options[i]);
		}

		if (iEndBasis != (options.size() - 1))
		{
			int iEndEcp;
			if (options[iEndBasis + 1] == "ActivateEcp")
			{
				for (size_t i = iEndBasis + 2; i < options.size(); i++)
				{
					if (options[i] == "EndOfEcp")
					{
						iEndEcp = i;
						break;
					}

					gamessEcpFiles.push_back(options[i]);
				}
			}
		}

		if ((gamessEcpFiles.size() != gamessAtomBasisFiles.size())
			&&
			(gamessEcpFiles.size()!= 0))
		{
			cout << "number of ecp and atomic bassis files don't match"
				<< endl << "check input" << endl;
			throw 0;
		}
		readGamessAuxFiles();
	}
}

void WriteQuantumInput::readGamessAuxFiles()
{
	bool ecp = gamessEcpFiles.size() != 0;
	int nAtoms = gamessAtomBasisFiles.size();
	atomBasis.resize(nAtoms);
	if (ecp)
		atomEcp.resize(nAtoms);

	for (int i = 0; i < nAtoms; i++)
	{
		ifstream readBasis_(("auxFiles/" + gamessAtomBasisFiles[i]).c_str());
		string auxline;
		while (getline(readBasis_, auxline))
		{
			stringstream line;
			string aux;
			line << auxline;
			line >> aux;
			if (aux == "end")
				break;

			atomBasis[i].push_back(auxline);
		}
		readBasis_.close();
		if (ecp)
		{
			ifstream readEcp_(("auxFiles/" + gamessEcpFiles[i]).c_str());
			string auxlineEcp;
			while (getline(readEcp_, auxlineEcp))
			{
				stringstream line;
				string aux;
				line << auxlineEcp;
				line >> aux;
				if (aux == "end")
					break;

				atomEcp[i].push_back(auxlineEcp);
			}
			readEcp_.close();
		}
	}
}

string WriteQuantumInput::createInput(
	vector<CoordXYZ> &coordinates,
	int indexAddedAtFinalOfInputToSaveItsName)
{
	// managing names
	stringstream convertInd;
	string inputName, indexString;
	convertInd << indexAddedAtFinalOfInputToSaveItsName;
	convertInd >> indexString;
	if (indexAddedAtFinalOfInputToSaveItsName == 0)
		indexString = "";

	inputName = projectName + indexString;
	try {
		if ((type == "mopac") || (type == "mopac2009"))
			buildMopacInput(coordinates, inputName);
		else if (type == "gamess")
			buildGamessInput(coordinates, inputName);
	}
	catch (...) {
		return "";
	}

	return inputName;
}

void WriteQuantumInput::buildGamessInput(vector<CoordXYZ> &coordinates, string inputName)
{
	inputName += ".inp";
	remove(inputName.c_str());
	ofstream gamessInput_(inputName.c_str());

	//writing header
	for (size_t i = 0; i < gamessHeader.size(); i++)
		gamessInput_ << gamessHeader[i] << endl;

	if (gamessAtomBasisFiles.size() != coordinates.size())
	{
		cout << "basis quantity and atoms quantity dont match" << endl;
		throw 0;
	}

	//writing coordinates and basis
	for (size_t i = 0; i < coordinates.size(); i++)
	{
		gamessInput_ << coordinates[i].atomlabel << "  "
			<< getChargeFromAtom(coordinates[i].atomlabel) << "  "
			<< coordinates[i].x << "  "
			<< coordinates[i].y << "  "
			<< coordinates[i].z << "  "
			<< endl;
		for (size_t j = 0; j < atomBasis[i].size(); j++)
			gamessInput_ << atomBasis[i][j] << endl;

		gamessInput_ << endl;
	}
	gamessInput_ << " $END" << endl;
	if (gamessEcpFiles.size() != 0)
	{
		gamessInput_ << " $ECP" << endl;
		for (size_t i = 0; i < atomEcp.size(); i++)
		{
			for (size_t j = 0; j < atomEcp[i].size(); j++)
				gamessInput_ << atomEcp[i][j] << endl;
		}

		gamessInput_ << " $END" << endl;
	}

	gamessInput_.close();
}

int WriteQuantumInput::getChargeFromAtom(string atomName)
{
	if ((atomName == "H")||(atomName == "h"))
		return 1;
	else if ((atomName == "He")||(atomName == "he")||(atomName == "HE"))
		return 2;
	else if ((atomName == "Li")||(atomName == "li")||(atomName == "LI"))
		return 3;
	else if ((atomName == "Be")||(atomName == "be")||(atomName == "BE"))
		return 4;
	else if ((atomName == "B")||(atomName == "b"))
		return 5;
	else if ((atomName == "C")||(atomName == "c"))
		return 6;
	else if ((atomName == "N")||(atomName == "n"))
		return 7;
	else if ((atomName == "O")||(atomName == "o"))
		return 8;
	else if ((atomName == "F")||(atomName == "f"))
		return 9;
	else if ((atomName == "Ne")||(atomName == "ne")||(atomName == "NE"))
		return 10;
	else if ((atomName == "Na")||(atomName == "na")||(atomName == "NA"))
		return 11;
	else if ((atomName == "Mg")||(atomName == "mg")||(atomName == "MG"))
		return 12;
	else if ((atomName == "Al")||(atomName == "al")||(atomName == "AL"))
		return 13;
	else if ((atomName == "Si")||(atomName == "si")||(atomName == "SI"))
		return 14;
	else if ((atomName == "P")||(atomName == "p"))
		return 15;
	else if ((atomName == "S")||(atomName == "s"))
		return 16;
	else if ((atomName == "Cl")||(atomName == "cl")||(atomName == "CL"))
		return 17;
	else if ((atomName == "Ar")||(atomName == "ar")||(atomName == "AR"))
		return 18;
	else if ((atomName == "K")||(atomName == "k"))
		return 19;
	else if ((atomName == "Ca")||(atomName == "ca")||(atomName == "CA"))
		return 20;
	else if ((atomName == "Sc")||(atomName == "sc")||(atomName == "SC"))
		return 21;
	else if ((atomName == "Ti")||(atomName == "ti")||(atomName == "TI"))
		return 22;
	else if ((atomName == "V")||(atomName == "v"))
		return 23;
	else if ((atomName == "Cr")||(atomName == "cr")||(atomName == "CR"))
		return 24;
	else if ((atomName == "Mn")||(atomName == "mn")||(atomName == "MN"))
		return 25;
	else if ((atomName == "Fe")||(atomName == "fe")||(atomName == "FE"))
		return 26;
	else if ((atomName == "Co")||(atomName == "co")||(atomName == "CO"))
		return 27;
	else if ((atomName == "Ni")||(atomName == "ni")||(atomName == "NI"))
		return 28;
	else if ((atomName == "Cu")||(atomName == "cu")||(atomName == "CU"))
		return 29;
	else if ((atomName == "Zn")||(atomName == "zn")||(atomName == "ZN"))
		return 30;
	else if ((atomName == "Ga")||(atomName == "ga")||(atomName == "GA"))
		return 31;
	else if ((atomName == "Ge")||(atomName == "ge")||(atomName == "GE"))
		return 32;
	else if ((atomName == "As")||(atomName == "as")||(atomName == "AS"))
		return 33;
	else if ((atomName == "Se")||(atomName == "se")||(atomName == "SE"))
		return 34;
	else if ((atomName == "Br")||(atomName == "br")||(atomName == "BR"))
		return 35;
	else if ((atomName == "Kr")||(atomName == "kr")||(atomName == "KR"))
		return 36;
	else if ((atomName == "Rb")||(atomName == "rb")||(atomName == "RB"))
		return 37;
	else if ((atomName == "Sr")||(atomName == "sr")||(atomName == "SR"))
		return 38;
	else if ((atomName == "Y")||(atomName == "y"))
		return 39;
	else if ((atomName == "Zr")||(atomName == "zr")||(atomName == "ZR"))
		return 40;
	else if ((atomName == "Nb")||(atomName == "nb")||(atomName == "NB"))
		return 41;
	else if ((atomName == "Mo")||(atomName == "mo")||(atomName == "MO"))
		return 42;
	else if ((atomName == "Tc")||(atomName == "tc")||(atomName == "TC"))
		return 43;
	else if ((atomName == "Ru")||(atomName == "ru")||(atomName == "RU"))
		return 44;
	else if ((atomName == "Rh")||(atomName == "rh")||(atomName == "RH"))
		return 45;
	else if ((atomName == "Pd")||(atomName == "pd")||(atomName == "PD"))
		return 46;
	else if ((atomName == "Ag")||(atomName == "ag")||(atomName == "AG"))
		return 47;
	else if ((atomName == "Cd")||(atomName == "cd")||(atomName == "CD"))
		return 48;
	else if ((atomName == "In")||(atomName == "in")||(atomName == "IN"))
		return 49;
	else if ((atomName == "Sn")||(atomName == "sn")||(atomName == "SN"))
		return 50;
	else if ((atomName == "Sb")||(atomName == "sb")||(atomName == "SB"))
		return 51;
	else if ((atomName == "Te")||(atomName == "te")||(atomName == "TE"))
		return 52;
	else if ((atomName == "I")||(atomName == "i"))
		return 53;
	else if ((atomName == "Xe")||(atomName == "xe")||(atomName == "XE"))
		return 54;
	else
	{
		cout << "ATOM NOT REGISTERED - CONTACT DEVELOPERS" << endl;
		throw 0;
	}
}










void WriteQuantumInput::buildMopacInput(
	vector<CoordXYZ> &allAtoms,
	string inputName)
{
	inputName += ".mop";
	remove(inputName.c_str());
	ofstream mopInput_(inputName.c_str());

	if (!(externalParams == ""))
		mopInput_ << "EXTERNAL=" + externalParams + ".inp  +" << endl;

	mopInput_ << mopacHeader << endl;
	mopInput_ << projectName << endl << endl;

	if (!(centralMetal == ""))
		mopInput_ << centralMetal << " 0.0   0   0.0   0   0.0   0 " << endl;

	for (size_t j = 0; j < allAtoms.size(); j++)
	{
		mopInput_ << allAtoms[j].atomlabel << "   "
			<< setprecision(16) << allAtoms[j].x << "  1  "
			<< setprecision(16) << allAtoms[j].y << "  1  "
			<< setprecision(16) << allAtoms[j].z << "  1  "
			<< endl;
	}
	mopInput_.close();
}


