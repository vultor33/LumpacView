#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "BuildComplex.h"
#include "RootMeanSquareDeviation.h"

using namespace std;

int main()
{
	/*
	--> TESTAR SE O RMS OVERLAY DA O MESMO RESULTADO DO AVOGADRO
	--> E TOCAR Ne

	Pegar o ligante original sem otimizar.
	Pegar o complexo de partida, la vai ter varios ligantes,
	  provavelmente diferentes - 1,2,3,4.
	Pegar o ligante que eu tenho interesse e criar um
	  ligante-lumpac-view dele.
	Eu vou construir uma geometria no input do lumpacviewinput
	  que se assemelha a que está lá - visualmente. Antes do
	  annealing.


	PASSO 2 - Colocar o metal no meio, preencher com agua e ou cloretos. 
	          - encontrar o posicionamento certo, entao vou precisar
			    recombinar o raio X até ter o menor rmsd possivel com o meu
				sobre os pontos na superficie.
	*/

	
	RootMeanSquareDeviation rmsd_;
	double rms = rmsd_.rmsOverlay("trip-antes-1.xyz", "trip-antes-2.xyz");
	cout << rms << endl;
	cin.get();
	// 0.053526610988479453
	
	/*
	int charge = 3;
	int coordination = 9;
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

	BuildComplex bc_;
	bc_.makeComplexOptimizingInMopac(ligandName, coordination, charge, options, mopacExecPath);
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









*/


