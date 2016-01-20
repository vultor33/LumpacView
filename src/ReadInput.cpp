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

ReadInput::ReadInput(){}

ReadInput::~ReadInput(){}

void ReadInput::rePrintInput()
{
	cout << "LUMPAC VIEW INPUT" << endl
		<< "metal name:     " << metalName << endl
		<< "metal params:   " << metalParams << endl
		<< "ligand files" << endl;
	for (int i = 0; i < (int)ligandFileName.size(); i++)
	{
		cout << ligandFileName[i] << endl;
	}
	cout << "END OF INPUT" << endl;
}


void ReadInput::readLumpacViewInput()
{
	ifstream input_(inputName.c_str());
	int nLigands;
	input_ >> metalName;
	input_ >> metalParams;
	input_ >> nLigands;
	ligandFileName.resize(nLigands);
	for (int i = 0; i < nLigands; i++)
	{
		input_ >> ligandFileName[i];
	}

	buildProjectName();
	readAllLigands();
}

void ReadInput::buildProjectName()
{
	projectName = metalName;
	for (int i = 0; i < (int)ligandFileName.size(); i++)
	{
		projectName += "-" + ligandFileName[i];
	}
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
	coord.resize(nAtoms);
	int i = 0;
	getline(mol_, auxline); //useless title
	while (getline(mol_, auxline))
	{
		if (auxline == "")
			break;

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
	//throw if x1 and x2 not found (close input)
	molecule.initializeLigand(coord); 
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


