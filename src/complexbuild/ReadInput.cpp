#include "ReadInput.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "MyExceptions.h"
#include "Ligand.h"
#include "Coordstructs.h"

using namespace std;

ReadInput::ReadInput()
{
	inputName = "LumpacViewInput.txt";
}

ReadInput::~ReadInput(){}

void ReadInput::readLumpacViewInput()
{
	// READING LUMPACVIEWINPUT
	ifstream input_(inputName.c_str());
	if (!input_.is_open()) { 
		MyExceptions mexept_(4);
		throw mexept_;
	}
	int nLigands;
	string metalName, metalParams;

	input_ >> metalName;
	input_ >> metalParams;
	input_ >> nLigands;
	ligandFileName.resize(nLigands);
	for (int i = 0; i < nLigands; i++)
	{
		input_ >> ligandFileName[i];
	}
	input_.close();

	options.resize(5);
	options[0] = "";
	options[1] = buildProjectName(metalName);
	options[2] = "";
	options[3] = metalParams;
	options[4] = metalName;

	readAllLigands();
}

void ReadInput::setLumpacViewInputWithoutFile(vector< string > & inputInformations)
{
	string metalName = inputInformations[0];
	string metalParams = inputInformations[1];
	for (size_t i = 0; i < (inputInformations.size() - 2); i++)
		ligandFileName.push_back(inputInformations[i+2]);

	options.resize(5);
	options[0] = "";
	options[1] = buildProjectName(metalName);
	options[2] = "";
	options[3] = metalParams;
	options[4] = metalName;

	readAllLigands();
}


void ReadInput::buildLumpacViewFromNames(std::vector< std::string > & fileNames)
{
	ligandFileName = fileNames;
	readAllLigands();
}


string ReadInput::buildProjectName(string metalName)
{
	string projName = metalName;
	for (int i = 0; i < (int)ligandFileName.size(); i++)
	{
		projName += "-" + ligandFileName[i];
	}
	return projName;
}

Ligand ReadInput::readConfigurations(string inputName)
{
	string nameinp = inputName + ".xyz";
	ifstream mol_(nameinp.c_str());
	if (!mol_.is_open()) 
	{ 
		MyExceptions mexept_(1);
		throw mexept_;
	}

	string auxline;
	int nAtoms;
	getline(mol_, auxline);
	stringstream line0;
	line0 << auxline;
	line0 >> nAtoms;

	vector<CoordXYZ> coord; 
	string titleInfo;
	coord.resize(nAtoms);

	int i = 0;
	getline(mol_, auxline);
	titleInfo = auxline;
	while (getline(mol_, auxline))
	{
		if (auxline == "")
			break;
		if (i == nAtoms) {
			MyExceptions mexept_("problem on ligand:  " + nameinp + " - check input");
			throw mexept_;
		}
		stringstream line;
		line << auxline;
		line >> coord[i].atomlabel
			>> coord[i].x
			>> coord[i].y
			>> coord[i].z;
		i++;
	}
	mol_.close();

	Ligand molecule;
	molecule.setLigandCoordinates(coord, titleInfo);
	return molecule;
}

void ReadInput::readAllLigands()
{
	int size = ligandFileName.size();
	allLigands.resize(size);
	for (int i = 0; i < size; i++)
	{
		allLigands[i] = readConfigurations(ligandFileName[i]);		
	}
}

void ReadInput::rePrintInput()
{
	cout << "LUMPAC VIEW INPUT" << endl
		<< "metal name:     " << options[4] << endl
		<< "metal params:   " << options[5] << endl
		<< "ligand files" << endl;
	for (int i = 0; i < (int)ligandFileName.size(); i++)
	{
		cout << ligandFileName[i] << endl;
	}
	cout << "END OF INPUT" << endl;
}



/*
if (!molecule.initializeLigand(coord)) {
MyExceptions mexept2_("problem on ligand:  " + nameinp + " - check input");
throw mexept2_;
}
*/
