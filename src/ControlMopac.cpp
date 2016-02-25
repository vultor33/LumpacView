#include "ControlMopac.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "Coordstructs.h"
#include "ReadInput.h"
#include "ReadMopac.h"
#include "Ligand.h"

using namespace std;

ControlMopac::ControlMopac(){}

ControlMopac::~ControlMopac(){}

bool ControlMopac::findLigandsConfiguration(ReadInput &ri_)
{
	optimization(ri_);

	return frequency > 0;
}


void ControlMopac::optimization(ReadInput &ri_)
{
	ReadMopac readmop_;
	Ligand auxLig;
	bool optimizeSucess;
	int optimizeCount = 0;

	Ligand optimizedMol_;
	do
	{
		pr_.buildMopacInput(ri_, all, "opt");

		system(("M2009_Ln_Orbitals.exe  " + ri_.projectName).c_str());

		Ligand auxOptimizedMol;
		optimizeSucess = readmop_.readOutput(ri_.projectName, auxOptimizedMol);

		if (optimizeSucess)
		{
			optimizedMol_ = auxOptimizedMol;
			ligands = all;
		}

		optimizeCount++;
		if (optimizeCount == 5)
		{
			cout << "otimization didn't succed - contact developers" << endl;\
			exit(1);
		}
	} while (!optimizeSucess);

	forceCalculation(optimizedMol_, ri_);
}


void ControlMopac::forceCalculation(Ligand &optimizedMol_, ReadInput &ri_)
{
	int k = 1;
	double centerX = optimizedMol_.coord[0].x;
	double centerY = optimizedMol_.coord[0].y;
	double centerZ = optimizedMol_.coord[0].z;
	for (int i = 0; i < (int)ligands.size(); i++) //n ligands
	{
		for (int j = 0; j < ((int)ligands[i].coord.size() - 2); j++) //n atoms on ligand i | removing x1 and x2
		{
			CoordXYZ auxCoord;
			auxCoord = optimizedMol_.coord[k];
			k++;
			ligands[i].coord[j].x = auxCoord.x - centerX;
			ligands[i].coord[j].y = auxCoord.y - centerY;
			ligands[i].coord[j].z = auxCoord.z - centerZ;
		}
	}

	pr_.buildMopacInput(ri_, ligands, "freq");

	system(("M2009_Ln_Orbitals.exe  " + ri_.projectName + "-freq").c_str());

	ReadMopac readmop_;
	if (!readmop_.readFrequency(ri_.projectName + "-freq"))
	{
		frequency = -1.0e0;
		cout << "frequency couldn't be evaluated" << endl;
	}
	else
		frequency = readmop_.getFirstFrequency();
}























void ControlMopac::buildMopacInput(
	vector<Ligand> &ligands, 
	string execOption)
{
	string header;
	string name;
	if (execOption == "opt"){
		name = projectName;
		header = mopacHeader;
	}
	else{
		name = projectName + "-freq";
		header = mopacFreq;
	}
	int allAtoms = 1;
	size_t size = ligands.size();
	for (size_t k = 0; k < size; k++)
		allAtoms += ligands[k].getNatoms();


	ofstream mopInput_(name.c_str());

	mopInput_ << "EXTERNAL=" + metalParams + ".inp  +" << endl;
	mopInput_ << header << endl;
	mopInput_ << projectName << endl << endl;
	mopInput_ << metalName << " 0.0   0   0.0   0   0.0   0 " << endl;
	vector<CoordXYZ> coord;
	for (size_t i = 0; i < ligands.size; i++)
	{
		coord = ligands[i].getAllAtoms();
		for (size_t j = 0; j < ligands[i].getNatoms(); j++)
		{
				mopInput_ << coord[j].atomlabel << "   "
					<< setprecision(16) << coord[j].x << "  1  "
					<< setprecision(16) << coord[j].y << "  1  "
					<< setprecision(16) << coord[j].z << "  1  "
					<< endl;
		}
	}
	mopInput_.close();
}
