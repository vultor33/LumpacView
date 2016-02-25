#include "ControlMopac.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

#include "Coordstructs.h"
#include "ReadMopac.h"
#include "Ligand.h"

using namespace std;

ControlMopac::ControlMopac(){}

ControlMopac::~ControlMopac(){}

bool ControlMopac::optimize(vector<Ligand> & ligands)
{
	vector<CoordXYZ> allAtoms;
	for (size_t i = 0; i < ligands.size(); i++)
	{
		vector<CoordXYZ> atomsLigand = ligands[i].getAllAtoms();
		allAtoms.insert(allAtoms.end(), atomsLigand.begin(), atomsLigand.end());
	}

	ReadMopac readmop_;

	buildMopacInput(allAtoms, "opt");

	system((mopacExecPath + projectName).c_str());

	vector<CoordXYZ> optimizedAtoms;

	optimizedAtoms = readmop_.readOutput(projectName);

	if (optimizedAtoms.size() == 0) return false;

	double frequency = forceCalculation(optimizedAtoms);

	return frequency > 0;

}

double ControlMopac::forceCalculation(vector<CoordXYZ> & optimizedAtoms)
{
	//metal always on 0 0 0.
	double centerX = optimizedAtoms[0].x;
	double centerY = optimizedAtoms[0].y;
	double centerZ = optimizedAtoms[0].z;
	vector<CoordXYZ> ligandAtoms(optimizedAtoms.size() - 1);
	for (size_t i = 0; i < optimizedAtoms.size(); i++)
	{
		ligandAtoms[i].x = optimizedAtoms[i + 1].x - centerX;
		ligandAtoms[i].y = optimizedAtoms[i + 1].y - centerY;
		ligandAtoms[i].z = optimizedAtoms[i + 1].z - centerZ;
	}

	buildMopacInput(ligandAtoms, "freq");

	system((mopacExecPath + projectName + "-freq").c_str());

	ReadMopac readmop_;

	double frequency = readmop_.readFrequency(projectName + "-freq");

	return frequency;
}

void ControlMopac::buildMopacInput(
	vector<CoordXYZ> &allAtoms, 
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
	ofstream mopInput_(name.c_str());

	mopInput_ << "EXTERNAL=" + metalParams + ".inp  +" << endl;
	mopInput_ << header << endl;
	mopInput_ << projectName << endl << endl;
	mopInput_ << metalName << " 0.0   0   0.0   0   0.0   0 " << endl;
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
