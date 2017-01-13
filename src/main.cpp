#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <time.h>
#include <math.h>

#include "BuildComplex.h"
#include "RootMeanSquareDeviation.h"
#include "BestPermutation.h"
#include "FindIsomers.h"
#include "ComplexCreator.h"
#include "AuxMath.h"
#include "DoAllPermutations.h"
#include "AllMolecularFormulas.h"
#include "CauchyIndex.h"

using namespace std;

void buildComplexWithALotOfIsomersAndDoWater(
	std::vector<int> permutation,
	std::string methodOptimize,
	std::string methodCosmo,
	std::string projectName
	)
{
	FindIsomers fd_;
	vector<string> options;
	string mopacExecPath;
	vector<CoordXYZ> allAtoms = fd_.buildComplexWithSelectedIsomer(
		permutation,
		projectName,
		methodOptimize,
		options,
		mopacExecPath);
	BuildComplex bc_;
	options[1] = projectName + "-water";
	options[2] = methodCosmo;
	options[3] = "";
	options[4] = "";
	bc_.runMopacAndPrint(options, mopacExecPath, allAtoms);
}

int main()
{
	CauchyIndex ci_;

	vector<CoordXYZ> mol(4);
	mol[0].x = 1.0e0;
	mol[0].y = 0.0e0;
	mol[0].z = 0.0e0;
	mol[1].x = 0.0e0;
	mol[1].y = 1.0e0;
	mol[1].z = 0.0e0;
	mol[2].x = -1.0e0;
	mol[2].y = 0.0e0;
	mol[2].z = 0.0e0;
	mol[3].x = 0.0e0;
	mol[3].y = -1.0e0;
	mol[3].z = 0.0e0;

	ci_.calculateAllIndexes();

	vector<int> permutation(5);
	vector<string> atoms(5);
	for (size_t i = 0; i < 5; i++)
	{
		permutation[i] = i;
		atoms[i] = "H ";
	}
	atoms[0] = "C ";
	atoms[1] = "C ";
	vector<int> bidentateAtomsChosen(2);
	bidentateAtomsChosen[0] = 0;
	bidentateAtomsChosen[1] = 1;

	ci_.printMolecule(permutation, atoms, bidentateAtomsChosen);




	// ATENCAO - CADASTRAR OS PONTOS HARDIN COM 12 CASAS

	/*
	RootMeanSquareDeviation rmsd_;
	double rms = rmsd_.rmsOverlay("zBUVXAR11.xyz", "zagua-rm1.xyz");
	cout << rms << endl;
	cin.get();
	*/

	/* MAKING COMPLEX
	int charge = 0;
	int coordination = 7;
	string ligandName = "JALNIU.xyz";
	string mopacExecPath = "M2009_Ln_Orbitals.exe";
	vector<string> options(5);
	options[0] = "mopac2009";
	options[1] = "JALNIU";
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

	/* 
	-> FIND ISOMERS - WORKING
	*/

//	FindIsomers findIso_;	
//	findIso_.start();

	/* CALCULANDO O COMPLEXO A PARTIR DOS CODIGOS 
	VERIFICAR - DOALLPERMUTATIONS
	vector<int> permutation(7);
	string projectName = "jalniu4";
	string methodOptimize = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	string methodCosmo = " EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
	permutation[0] = 0;	permutation[1] = 1;	permutation[2] = 2;	permutation[3] = 3;	permutation[4] = 4;	permutation[5] = 5;	permutation[6] = 6;
	buildComplexWithALotOfIsomersAndDoWater(permutation, methodOptimize, methodCosmo, projectName);
	*/


//  COMPLEX WITH SELECTED ISOMER

/*
	FindIsomers fd_;
	vector<string> options;
	vector<int> permutation(7);
	string projectName;
	string method;
	string mopacExecPath;
	permutation[0] = 0;
	permutation[1] = 1;
	permutation[2] = 2;
	permutation[3] = 3;
	permutation[4] = 4;
	permutation[5] = 5;
	permutation[6] = 6;
	projectName = "jalniu3";
	method = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	vector<CoordXYZ> res1 = fd_.buildComplexWithSelectedIsomer(
		permutation,
		projectName,
		method,
		options,
		mopacExecPath);
	BuildComplex bc_;
	options[1] = projectName + "-water";
	options[2] = " EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
	options[3] = "";
	options[4] = "";
	bc_.runMopacAndPrint(options, mopacExecPath, res1);
*/

/*
	DoAllPermutations doall_;
	string methodOptimize = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	string methodCosmo = " RM1 EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
	int nPermutations = 7;
	string projectName = "jalniu-RM1-orbitais";
	string filePermutations = "BBm1m1m1.xyz";
	doall_.calculateAll(
		methodOptimize,
		methodCosmo,
		projectName,
		filePermutations,
		nPermutations);
*/

/*
	ifstream all_("nomes.txt");
	string auxline;
	while (getline(all_, auxline))
	{
		string name;
		stringstream line;
		line << auxline;
		line >> name;
		if (name == "end")
			break;
		else
		{
			DoAllPermutations doall;
			doall.analysis("!JALNIU-CRISTAL.xyz", name.c_str());
		}
	}
*/

//	RootMeanSquareDeviation rmsd_;


//	AllMolecularFormulas comb_;
//	comb_.doAllCombinations(1);
	//comb_.findAllIsomersOnCombinations("allCombinations.txt");
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



/*
DoAllPermutations doall_;
string methodOptimize = " AM1 SPARKLE BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo = " AM1 SPARKLE EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
int nPermutations = 7;
string projectName = "jalniu-AM1";
string filePermutations = "BBm1m1m1.xyz";
doall_.calculateAll(
methodOptimize,
methodCosmo,
projectName,
filePermutations,
nPermutations);

DoAllPermutations doall2_;
string methodOptimize2 = " RM1 SPARKLE BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo2 = " RM1 SPARKLE EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
string projectName2 = "jalniu-RM1";
doall2_.calculateAll(
methodOptimize2,
methodCosmo2,
projectName2,
filePermutations,
nPermutations);

DoAllPermutations doall3_;
string methodOptimize3 = " PM3 SPARKLE BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo3 = " PM3 SPARKLE EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
string projectName3 = "jalniu-PM3";
doall3_.calculateAll(
methodOptimize3,
methodCosmo3,
projectName3,
filePermutations,
nPermutations);

DoAllPermutations doall4_;
string methodOptimize4 = " PM7 SPARKLE BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo4 = " PM7 SPARKLE EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
string projectName4 = "jalniu-PM7";
doall_.calculateAll(
methodOptimize4,
methodCosmo4,
projectName4,
filePermutations,
nPermutations);

DoAllPermutations doall5_;
string methodOptimize5 = " PM6 SPARKLE BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo5 = " PM6 SPARKLE EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
string projectName5 = "jalniu-PM6";
doall_.calculateAll(
methodOptimize5,
methodCosmo5,
projectName5,
filePermutations,
nPermutations);

DoAllPermutations doall6_;
string methodOptimize6 = " PM6-D3H4 SPARKLE BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo6 = " PM6-D3H4 SPARKLE EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
string projectName6 = "jalniu-PM6-D3H4";
doall_.calculateAll(
methodOptimize6,
methodCosmo6,
projectName6,
filePermutations,
nPermutations);
*/


/* CODIGO PARA VERIFICAR O TEMPO DE UMA ROTINA

clock_t begin = clock();

rotina AQUI

clock_t end = clock();
double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

cout << "demorou:  " << elapsed_secs << endl;
*/




