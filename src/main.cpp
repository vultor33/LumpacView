#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "BuildComplex.h"

using namespace std;

int main()
{
	/*
	PASSO 1 - construir todos os ligantes no formato que está aqui.


	PASSO 2 - Colocar o metal no meio, preencher com agua e ou cloretos. 
	          - encontrar o posicionamento certo, entao vou precisar
			    recombinar o raio X até ter o menor rmsd possivel com o meu
				sobre os pontos na superficie.
	*/

	int charge = 3;
	int coordination = 9;
	string ligandName = "h2o-coord-9.xyz";
	string mopacExecPath = "M2009_Ln_Orbitals.exe";
	vector<string> options(5);
	options[0] = "mopac2009";
	options[1] = "h2o-teste-mop";
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

	return 0;
}




