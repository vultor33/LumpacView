#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include "BuildComplex.h"
#include "RootMeanSquareDeviation.h"

using namespace std;

void printCoordXYZ(vector<CoordXYZ> & allAtoms, string fName)
{
	ofstream pr_(fName.c_str());
	pr_ << allAtoms.size() << endl << "useless" << endl;
	for (size_t i = 0; i < allAtoms.size(); i++)
	{
		pr_ << allAtoms[i].atomlabel << "  "
			<< allAtoms[i].x << "  "
			<< allAtoms[i].y << "  "
			<< allAtoms[i].z << endl;
	}

	pr_.close();
}


vector<CoordXYZ> ligandToCoordXYZ(vector<Ligand> & allLigands)
{
	vector<CoordXYZ> newAllAtoms(1);
	newAllAtoms[0].atomlabel = "Eu";
	newAllAtoms[0].x = 0.0e0;
	newAllAtoms[0].y = 0.0e0;
	newAllAtoms[0].z = 0.0e0;
	for (size_t i = 0; i < allLigands.size(); i++)
	{
		vector<CoordXYZ> ligandAdd = allLigands[i].getAllAtoms();
		newAllAtoms.insert(
			newAllAtoms.end(),
			ligandAdd.begin(),
			ligandAdd.end());
	}
	return newAllAtoms;
}



unsigned int factorial(unsigned int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}

vector< vector<int> > allFactorialPermutations(const int nMax)
{
	int * myints;
	myints = new int[nMax];
	for (int i = 0; i < nMax; i++)
		myints[i] = i+1;

	std::sort(myints, myints + nMax);

	vector< vector<int> > permutations(factorial(nMax));
	for (int i = 0; i < permutations.size(); i++)
		permutations[i].resize(nMax);

	int k = 0;
	do
	{		
		permutations[k].resize(nMax);
		for (int i = 0; i < nMax; i++)
		{
			permutations[k][i] = myints[i];
			cout << myints[i] << "  ";
		}
		cout << endl;
		k++;
	} while (std::next_permutation(myints, myints + nMax));

//	std::cout << "After loop: " << myints[0] << ' ' << myints[1] << ' ' << myints[2] << ' ' << myints[3] << '\n';

	delete[] myints;

	return permutations;
}


















int main()
{
	/*

	SUGESTAO ---> A OTIMIZACAO DEVERIA SER FEITA COM BFGS	             

	*/

	/*
	RootMeanSquareDeviation rmsd_;
	double rms = rmsd_.rmsOverlay("zBUVXAR11.xyz", "zagua-rm1.xyz");
	cout << rms << endl;
	cin.get();
	*/

	/* MAKING COMPLEX
	int charge = 1;
	int coordination = 8;
	string ligandName = "DUCNAQ-OONO.xyz";
	string mopacExecPath = "M2009_Ln_Orbitals.exe";
	vector<string> options(5);
	options[0] = "mopac2009";
	options[1] = "DUCNAQ-OONO";
	options[2] = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n"; //freq = AUX THERMO FORCE
																		//adding charge
	string chargeString;
	stringstream convert;
	convert << charge;
	convert >> chargeString;
	chargeString = " CHARGE=" + chargeString + " ";
	// charge computed
	options[2] += " NOLOG GEO-OK SCFCRT=1.D-10" + chargeString;
	options[3] = "Eu_spk";
	options[4] = "Eu";
	BuildComplex bc_;
	bc_.makeComplexOptimizingInMopac(ligandName, coordination, charge, options, mopacExecPath);
	*/

	/* RUNNING THROUGH LUMPAC VIEW INPUT
	BuildComplex bc_;
	vector<CoordXYZ> vec = bc_.build();
	*/

	// partindo do lumpac view input temos aqui o primeiro assemble
	RootMeanSquareDeviation rmsd_;
	vector<CoordXYZ> molCrystal = rmsd_.readCoord("DUCNAQ.xyz");

	BuildComplex bc_;
	vector<Ligand> allAtoms = bc_.assembleComplexWithoutSA();

	printCoordXYZ(ligandToCoordXYZ(allAtoms),"assembleLigands.xyz");

	printCoordXYZ(molCrystal, "znormal-antes-cristal.xyz");

	vector<CoordXYZ> mol2 = ligandToCoordXYZ(allAtoms);

	double rms = rmsd_.rmsOverlay(molCrystal,mol2);

	molCrystal.insert(molCrystal.end(), mol2.begin(), mol2.end());

	printCoordXYZ(molCrystal, "znormal-depois-super.xyz");


	vector< vector<int> > allPerm = allFactorialPermutations(4);



	cin.get();


	/*
	- na verdade eu preciso dos ligantes nesse ponto aqui.
	- vou ordena-los de varias formas e comparar com o kabsch algorithm
	*/

	return 0;
}

/*
CRIACAO DE LIGANTES

int charge = 3;
int coordination = 9;
string ligandName = "BUVXAR11.xyz";
string mopacExecPath = "M2009_Ln_Orbitals.exe";
vector<string> options(5);
options[0] = "mopac2009";
options[1] = "buvxar11";
options[2] = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n"; //freq = AUX THERMO FORCE
//adding charge
string chargeString;
stringstream convert;
convert << charge;
convert >> chargeString;
chargeString = " CHARGE=" + chargeString + " ";
// charge computed
options[2] += " NOLOG GEO-OK SCFCRT=1.D-10" + chargeString;
options[3] = "Eu_spk";
options[4] = "Eu";
BuildComplex bc_;
bc_.makeComplexOptimizingInMopac(ligandName, coordination, charge, options, mopacExecPath);

DUCNAQ
int charge = 1;
int coordination = 8;
string ligandName = "DUCNAQ-ligand.xyz";
string mopacExecPath = "M2009_Ln_Orbitals.exe";
vector<string> options(5);
options[0] = "mopac2009";
options[1] = "DUCNAQ-ligand";
options[2] = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n"; //freq = AUX THERMO FORCE
//adding charge
string chargeString;
stringstream convert;
convert << charge;
convert >> chargeString;
chargeString = " CHARGE=" + chargeString + " ";
// charge computed
options[2] += " NOLOG GEO-OK SCFCRT=1.D-10" + chargeString;
options[3] = "Eu_spk";
options[4] = "Eu";

TESTE DO KABSCH RMSD
RootMeanSquareDeviation rmsd_;
double rms = rmsd_.rmsOverlay("trip-antes-1.xyz", "trip-antes-2.xyz");
cout << rms << endl;
cin.get();
*/


