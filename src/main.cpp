#include <utility>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <time.h>
#include <math.h>
//#include <unistd.h>
#include <sys/stat.h>

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

inline bool exist_file (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


int main(int argc, char *argv[])
{
	clock_t begin = clock();

	CauchyIndex ci_(6);
	vector<int> permut(6);
	for (int i = 0; i < 6; i++)
		permut[i] = i;
	ci_.enantiomersOrderingBlock(2, 30, "skeleton-6.txt", permut);
	return 0;


	vector<string> atom;
	vector<int> bidentate;
	ci_.enantiomersOrdering();
	return 0;

	// TENHO QUE CRIAR UMA PARADA QUE SEPARE OS ISOMEROS NA GRAFENO.
	// acho que e correr todos com todos mesmo, fazer o q.
	// de cima pra baixo com numeros lidos atraves do modulo.

	// colocar todos para trabalhar juntos, pega o primeiro
	// divide o resto nos processadores (80 mil em cada)
	// fica tipo - do 400 ao 700.
	// se nao tiver nada eu descarto. (quando tiver eu bato la)







	// ATENCAAAAAOOOOOOOOOOOOOOOOOOOO
	// a comparacao de isomeros do 12 esta restrita

/*	
	string flagsFile = "results-m01m02m03m04m04m04";
	FindIsomers fd_;
	vector<int> permutation(8);
	permutation[0] = 0;
	permutation[1] = 1;
	permutation[2] = 2;
	permutation[3] = 3;
	permutation[4] = 4;
	permutation[5] = 6;
	permutation[6] = 5;
	permutation[7] = 7;
	fd_.printSelectedIsomerCauchy(permutation, flagsFile);
	return 0;
	*/

// !!! IMPORTANTE - mudar bool ComplexCreator::start(vector<int> & ligandsPermutation)

// ele troca as coordenadas




/*    	 EXEMPLO PARA RODAR
	./lumpacview.exe generateCompositionFiles 6 30 12 3 m01m02m03m04m04m04 skeleton-6 pc
	./lumpacview.exe cleanBlocks m01m02m03m04m04m04 6 3 pc
	cat composition---atomTypes composition-independent* > results-composition

	criar uma function que retorna o working directory
*/

	stringstream convert0;
	convert0 << argv[1];
	string execType = convert0.str();
	if(execType == "wholeBlockDeletion")
	{
        int systemSize, kRotateInit, kRotateEnd, lDeleteInit, lDeleteEnd;
        stringstream cGen;
        cGen << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5] << "  " << argv[6];
        cGen >> systemSize >> kRotateInit >> kRotateEnd >> lDeleteInit >> lDeleteEnd;
		CauchyIndex ci_(systemSize);
		if(systemSize == 12)
			ci_.doBlockRAMDeletion12(kRotateInit, kRotateEnd, lDeleteInit, lDeleteEnd);
		else
			ci_.doBlockRAMDeletion(kRotateInit, kRotateEnd, lDeleteInit, lDeleteEnd);
	}
	else if(execType == "generateFiles")
	{
        int systemSize, nProc;
		string machineType;
        stringstream cGen;
        cGen << argv[2] << "  " << argv[3] << "  " << argv[4];
        cGen >> systemSize >> nProc >> machineType;
		CauchyIndex ci_(systemSize);
		ci_.generateSlurmFilesToDeletion(systemSize, nProc, machineType);
	}
	else if(execType == "runall")
	{
	    system("pwd > workingDir");
       	ifstream wDir_("workingDir");
	    string workingDirectory;
        wDir_ >> workingDirectory;
	    wDir_.close();
		remove("workingDir");
        int systemSize, blockInit, nProc;
		string machineType;
        stringstream cGen;
        cGen << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5];
        cGen >> systemSize >> blockInit >> nProc >> machineType;
		CauchyIndex ci_(systemSize);
		if(blockInit <= nProc)
			ci_.runall(blockInit, nProc,machineType,workingDirectory);
	}
	else if(execType == "blockGeneration")
	{
		int systemSize;
	        string blockName;
		string composition;
        	stringstream cGen;
        	cGen << argv[2] << "  " << argv[3] << "  " << argv[4];
        	cGen >> composition >> systemSize >> blockName;
        	CauchyIndex ci_(systemSize);
		if(composition == "none")
			ci_.generateAllIndependentIsomersRuntimeRotationsAndReadBlock(blockName);
		else
			ci_.generateAllIndependentIsomersWithFlag(blockName, composition + "---atomTypes.txt", composition);
	}
	else if(execType == "cleanBlocks")
	{
		system("pwd > workingDir");
		ifstream wDir_("workingDir");
		string workingDirectory;
		wDir_ >> workingDirectory;
		wDir_.close();
		remove("workingDir");
        	int systemSize, nMaxBlock;
		string composition;
		string machineType;
        	stringstream cGen;
        	cGen << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5];
        	cGen >> composition >> systemSize >> nMaxBlock >> machineType;
		CauchyIndex ci_(systemSize);
		ci_.cleanBlocksAndGenerateIsomers(nMaxBlock,systemSize,composition,workingDirectory,machineType);

	}
	else if(execType == "generateCompositionFiles")
	{
		system("pwd > workingDir");
                ifstream wDir_("workingDir");
                string workingDirectory;
                wDir_ >> workingDirectory;
                wDir_.close();
                remove("workingDir");

		int systemSize, total, bigBlockSize, smallBlockSize;
		string composition, rawIsomersFile, machineType;

		stringstream cGen;
		cGen << argv[2] << "  " 
			<< argv[3] << "  " 
			<< argv[4] << "  " 
			<< argv[5] << "  " 
			<< argv[6] << "  " 
			<< argv[7] << "  " 
			<< argv[8];
		
		cGen >> systemSize
			>> total 
			>> bigBlockSize
			>> smallBlockSize
			>> composition
			>> rawIsomersFile
			>> machineType;
		
		CauchyIndex ci_(systemSize);
		ci_.generateAtomTypesAndBidentateChosenFile(composition);
		
		string compositionFile = workingDirectory + "/" + composition + "---atomTypes.txt";

		ci_.generateSlurmFilesToDeletionFlags(
			systemSize,
			total,
			bigBlockSize,	
			smallBlockSize,
			compositionFile,
			workingDirectory + "/" + rawIsomersFile,
			workingDirectory,
			machineType);


	}	
	else if (execType == "compositionBlockDeletion")
	{
		int systemSize, kRotateInit, kRotateEnd, lDeleteInit, lDeleteEnd;
		stringstream cGen;
		string rawIsomers, compositionFile;
		cGen << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5] << "  " << argv[6] << "  " << argv[7] << "  " << argv[8];
		cGen >> systemSize >> rawIsomers >> compositionFile >> kRotateInit >> kRotateEnd >> lDeleteInit >> lDeleteEnd;
		CauchyIndex ci_(systemSize);
		ci_.doBlockDeletionFlags(
			rawIsomers,
			compositionFile,
			kRotateInit,
			kRotateEnd,
			lDeleteInit,
			lDeleteEnd);
	}
	else
		cout << "execType not found" << endl;


	/* FUNCIONANDO PARA DELETAR OS CARAS
        int systemSize, kRotateInit, kRotateEnd, lDeleteInit, lDeleteEnd;
        stringstream cGen;
        cGen << argv[1] << "  " << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5];
        cGen >> systemSize >> kRotateInit >> kRotateEnd >> lDeleteInit >> lDeleteEnd;
	CauchyIndex ci_(systemSize);
	ci_.doBlockRAMDeletion(kRotateInit, kRotateEnd, lDeleteInit, lDeleteEnd);
	*/

	/* TESTANDO PARA DELETAR O CARA
	int systemSize = 9;
	vector<int> permutation(systemSize);
        for (size_t i = 0; i < systemSize; i++)
                permutation[i] = i;
	vector< vector<int> > allPermutations;
        do
        {
		allPermutations.push_back(permutation);

	} while (std::next_permutation(permutation.begin(), permutation.end()));

	vector<int> apaga(systemSize);
	for(size_t i = 0; i < apaga.size(); i++)
	{
		if(i == apaga.size() - 1)
			apaga[i] = 0;
		else
			apaga[i] = i + 1;
	}

	vector< vector<int> > auxPerm;
	auxPerm = allPermutations;
	for(size_t i = 0; i < 100; i++)
	{
		auxPerm = allPermutations;


		for(size_t j = 0; j < auxPerm.size(); j++)
		{
			if(apaga == auxPerm[j])
			{
				vector< vector<int> >::iterator it = auxPerm.begin() + j;
				rotate(it, it+1,auxPerm.end());
				auxPerm.pop_back();
				break;
			}

		}


		//auxPerm = allPermutations;
		//auxPerm.erase(std::remove(auxPerm.begin(), auxPerm.end(), apaga), auxPerm.end());
		//inves de usar begin e end eu posso partir isso em 100 usar begin + 100, begin +200 e asism vai.
		//e soltar em todos os processafores
	}
	*/

	/*
	for(size_t i = 0; i < auxPerm.size(); i++)
	{
		for(size_t j = 0; j < systemSize; j++)
		{
			cout << auxPerm[i][j] << "  ";
		}
		cout << endl;
	}
	*/




	/*
	string execType;
	int systemSize, kInit, kEnd;
	stringstream cGen;
	cGen << argv[1] << "  " << argv[2] << "  " << argv[3] << "  " << argv[4];
	cGen >> execType >> systemSize >> kInit >> kEnd;
	CauchyIndex ci_(systemSize);
	if(execType == "genBlock")
		ci_.generateBlockFiles(systemSize,kInit,kEnd);
	else if(execType == "deletion")
		ci_.doBlockDeletion(kInit,kEnd);
	*/

	//ci_.doBlockDeletion(1,20);

	//ci_.generateBlockFiles(6,700,720);

	//string blockName = argv[1];

	//ci_.generateAllIndependentIsomers12(blockName);

	//ci_.generateAllIndependentIsomersRuntimeRotationsAndReadBlock(blockName);

	//ci_.generateAllIndependentIsomersRuntimeRotations();

	/*
	vector<string> allBlockNames(12);
	allBlockNames[0] = "isomersOf-block-10---1.txt";
	allBlockNames[1] = "isomersOf-block-10---7.txt";
	allBlockNames[2] = "isomersOf-block-10---2.txt";
	allBlockNames[3] = "isomersOf-block-10---8.txt";
	allBlockNames[4] = "isomersOf-block-10---3.txt";
	allBlockNames[5] = "isomersOf-block-10---9.txt";
	allBlockNames[6] = "isomersOf-block-10---4.txt";
	allBlockNames[7] = "isomersOf-block-10---10.txt";
	allBlockNames[8] = "isomersOf-block-10---5.txt";
	allBlockNames[9] = "isomersOf-block-10---11.txt";
	allBlockNames[10] = "isomersOf-block-10---6.txt";
	allBlockNames[11] = "isomersOf-block-10---12.txt";
	ci_.mergeBlocks(allBlockNames,12);
	*/


//	vector<string> atoms;
//	vector<int> bidentaChos;
//	ci_.generateAllIndependentIsomersRuntimeRotations();
	//ci_.rotationTest(atoms, bidentaChos);

/*
	vector<int> bidentateAtomsChosen;
	vector<int> atomTypes(6);
	vector<int> permutationIsomer1(6);
	vector<int> permutationIsomer2(6);
	for (size_t i = 0; i < 6; i++)
	{
		atomTypes[i] = 0;
		permutationIsomer1[i] = i;
		permutationIsomer2[i] = i;
	}
	permutationIsomer2[1] = 2;
	permutationIsomer2[2] = 3;
	permutationIsomer2[3] = 4;
	permutationIsomer2[4] = 1;
	//	atomTypes[3] = 0;
	bidentateAtomsChosen.resize(6);
	bidentateAtomsChosen[0] = 0;
	bidentateAtomsChosen[1] = 1;
	bidentateAtomsChosen[2] = 2;
	bidentateAtomsChosen[3] = 3;
	bidentateAtomsChosen[4] = 4;
	bidentateAtomsChosen[5] = 5;
	int compare = ci_.compareTwoIsomers(
		atomTypes,
		bidentateAtomsChosen,
		permutationIsomer1,
		permutationIsomer2);

*/

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "demorou:  " << elapsed_secs << "  segundos" << endl;

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

	/* CHECAR LumpacViewInput.txt para informacoes do input
	clock_t begin = clock();
	FindIsomers findIso_;
	findIso_.start();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "demorou:  " << elapsed_secs << "  segundos" << endl;
	*/

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




