#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <time.h>
#include <math.h>

#include "BuildComplex.h"
#include "RootMeanSquareDeviation.h"
#include "BestPermutation.h"
#include "FindIsomers.h"


using namespace std;


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

	// RUNNING THROUGH LUMPAC VIEW INPUT
	//BuildComplex bc_;
	//vector<CoordXYZ> vec = bc_.build();

	/* BUXVAR11
	vector<string> ligandNames(9);
	ligandNames[0] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[1] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[2] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[3] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[4] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[5] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[6] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[7] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[8] = "Lumpac-View-Ligand-BUVXAR11";
	BestPermutation bp_(ligandNames, "BUVXAR11.xyz");
	*/

	
	/* DUCNAQ WORKING
	vector<string> ligandNames(5);
	ligandNames[0] = "Lumpac-View-Ligand-DUCNAQ-ligand-FALTA-OTIMIZAR";
	ligandNames[1] = "Lumpac-View-Ligand-DUCNAQ-OONO";
	ligandNames[2] = "Lumpac-View-Ligand-BUVXAR11";
	ligandNames[3] = "Lumpac-View-Ligand-DUCNAQ-OONO";
	ligandNames[4] = "Lumpac-View-Ligand-BUVXAR11";
	BestPermutation bp_(ligandNames, "DUCNAQ.xyz");
	bp_.findBestPermutation();
	*/

	/* VIGPAC - WORKING
	vector<string> ligandNames(7);
	ligandNames[0] = "ROCTAF-CL";
	ligandNames[1] = "ROCTAF-CL";
	ligandNames[2] = "ROCTAF-CL";
	ligandNames[3] = "ROCTAF-ring";
	ligandNames[4] = "ROCTAF-ring";
	ligandNames[5] = "ROCTAF-ring";
	ligandNames[6] = "ROCTAF-ring";
	BestPermutation bp_(ligandNames, "VIGPAC-teste.xyz");
	bp_.findBestPermutation();
	*/

	/* ROCTAF
	vector<string> ligandNames(5);
	ligandNames[0] = "ROCTAF-alif";
	ligandNames[1] = "ROCTAF-ring";
	ligandNames[2] = "ROCTAF-CL";
	ligandNames[3] = "ROCTAF-CL";
	ligandNames[4] = "ROCTAF-CL";
	BestPermutation bp_(ligandNames, "roctaf-complex.xyz");
	bp_.findBestPermutation();
	*/


	/* SERHED  
	vector<string> ligandNames(3);
	ligandNames[0] = "SERHED-ligand-ciclo-5";
	ligandNames[1] = "SERHED-ligand-ciclo-5";
	ligandNames[2] = "SERHED-mono";
	*/
	
	
	//BestPermutation bp_("serhed-complex.xyz");
	//bp_.findBestPermutation();
	

	/* SOPFUY 
	vector<string> ligandNames(3);
	ligandNames[0] = "SOPFUY-ligand";
	ligandNames[1] = "SOPFUY-ligand";
	ligandNames[2] = "SOPFUY-ligand";
	BestPermutation bp_(ligandNames, "sopfuy-complex.xyz");
	bp_.findBestPermutation();
	*/

	/* SUXXIS 
	vector<string> ligandNames(6);
	ligandNames[0] = "SUXXIS-ligand";
	ligandNames[1] = "SUXXIS-ligand";
	ligandNames[2] = "SUXXIS-ligand";
	ligandNames[3] = "SUXXIS-ligand";
	ligandNames[4] = "SUXXIS-ligand";
	ligandNames[5] = "SUXXIS-ligand";
	BestPermutation bp_(ligandNames, "suxxis-complex.xyz");
	bp_.findBestPermutation();  */

	//  Start Timers
	//double wall0 = clock();

	FindIsomers findIso_;
	findIso_.start();

	//  Stop timers
	//double wall1 = clock();

	//cout << "Wall Time = " << wall1 - wall0 << endl;

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


