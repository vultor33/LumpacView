//#define UNIX

#include <utility>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <time.h>d
#include <math.h>
#include <sys/stat.h>
#include <iomanip>
#ifdef UNIX
	#include <unistd.h>
#endif

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

void waitSlurmFinish(int usedProc);

void buildCsvFile(int size, string skeletonName);

void changeNameOfFiles(string name);

string sizeToGeometryCode(int size);

int main(int argc, char *argv[])
{
	/*
	string line;
	string responseName = "response-combinations4.txt";
	ifstream resp_(responseName.c_str());
	while (!resp_.eof())
	{
		getline(resp_, line);
		if (line == "")
			break;
		stringstream convertLine;
		convertLine << line;
		string code;
		convertLine >> code;
		CauchyIndex ci123_(4);
		ci123_.generateAtomTypesAndBidentateChosenFile(code);
		ci123_.generateAllIndependentIsomersWithFlagEnantiomers("enatiomers-4.log", code + "---atomTypes.txt", code);
	}
	string code = "m01m02m03m04m05m06";
	CauchyIndex ci123_(6);
	ci123_.generateAtomTypesAndBidentateChosenFile(code);
	ci123_.generateAllIndependentIsomersWithFlagEnantiomers("enantiomers-6.log", code + "---atomTypes.txt", code);
	exit(0);
	*/

	CauchyIndex ci234(6);
	vector<int> permutation(6);
	vector<int> atomTypes(6);
	for (int i = 0; i < 6; i++)
	{
		permutation[i] = i;
		atomTypes[i] = i;
	}
	atomTypes[1] = 0;
	atomTypes[2] = 0;
	atomTypes[3] = 0;
	permutation[4] = 1;
	permutation[1] = 4;
	vector<int> bidChosen;
	ci234.indetifyIsomer(permutation, atomTypes, bidChosen);

	string responseName;
	cout << "type line: " << endl;
	//cin >> responseName;
	responseName = "response-combinations8.txt";
	changeNameOfFiles(responseName);
	return 0;

	clock_t begin = clock();

	//CauchyIndex ci123_(10);
	//ci123_.temporario();
	//return 0;
	/*
	string rawIsomers2 = "enantiomers-6";
	string composition3 = "m01m02m03m04m05m05";
	CauchyIndex ci14_(6);
	int kInit3 = 1;
	int kFinal3 = 10;
	int lDeleteI = 11;
	int lDeleteF = 30;
	ci14_.generateAtomTypesAndBidentateChosenFile(composition3);
	ci14_.doBlockDeletionFlagsEnantiomers(
			rawIsomers2,
                        composition3 + "---atomTypes.txt" ,
                        kInit3,
                        kFinal3,
                        lDeleteI,
                        lDeleteF);
	return 0;

	int systemSize = 7;
	string composition = "m01m02m02B01C01";
	string blockName = "enantiomers-7.log";
	CauchyIndex ci13_(systemSize);
	ci13_.printAllMoleculesFromFile(composition);
	return 0;
	ci13_.generateAtomTypesAndBidentateChosenFile(composition);
	ci13_.generateAllIndependentIsomersWithFlagEnantiomers(blockName, composition + "---atomTypes.txt", composition);
	//blockName = "skeleton-6.txt";
	//ci13_.generateAllIndependentIsomersWithFlag(blockName, composition + "---atomTypes.txt", composition);
	return 0;
	*/
	/*
	CauchyIndex ci5_(6);
	ci5_.generateAllIndependentIsomers();
	ci5_.createEnantiomersFiles(
		4,
		30,
		"skeleton-6",
		"",
		"pc");
	return 0;

	CauchyIndex ci_(6);
	int iPer = 7;
	ci_.enantiomersOrderingBlock(iPer, 30, "skeleton-6.txt");
	return 0;
*/
	// agora tenho que fazer aquele negocio ---
	// os numeros precisam variar de 1 - 31 - 61
	// e assim vai. cada processador vai fazer uma serie dessas.
	// definida a serie e so chamar
	// ci_.enantiomersOrderingBlock(i + 1, iFinal, "skeleton-6.txt", permut(i));
	// que resolve

/*
	vector<string> atom;
	vector<int> bidentate;
	ci_.enantiomersOrdering();
	return 0;
*/
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




/*  EXEMPLO PARA RODAR
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
			ci_.generateAllIndependentIsomersWithFlagEnantiomers(blockName, composition + "---atomTypes.txt", composition);
			//ci_.generateAllIndependentIsomersWithFlag(blockName, composition + "---atomTypes.txt", composition);
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
		//ci_.doBlockDeletionFlags(
		ci_.doBlockDeletionFlagsEnantiomers(
			rawIsomers,
			compositionFile,
			kRotateInit,
			kRotateEnd,
			lDeleteInit,
			lDeleteEnd);
	}
	else if (execType == "generateEnantiomersFiles")
	{
                system("pwd > workingDir");
                ifstream wDir_("workingDir");
                string workingDirectory;
                wDir_ >> workingDirectory;
                wDir_.close();
                remove("workingDir");
                int systemSize, nProc, iMax;
		string skeletonFile;
                string machineType;
                stringstream cGen;
                cGen << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5] << "  " << argv[6];
                cGen >> systemSize >> nProc >> iMax >> skeletonFile >> machineType;
                CauchyIndex ci_(systemSize);
        	ci_.createEnantiomersFiles(
                	nProc,
                	iMax,
                	skeletonFile,
                	workingDirectory,
                	machineType);
	}
	else if (execType == "enantiomerDeletion")
	{
		int systemSize, iPer, iMax, nProc;
		stringstream cGen;
		string skeletonFile;
		cGen << argv[2] << "  " << argv[3] << "  " << argv[4] << "  " << argv[5] << "  " << argv[6];
		cGen >> systemSize >> iPer >> iMax >> nProc >> skeletonFile;
		CauchyIndex ci_(systemSize);
		ci_.enantiomersOrderingBlock(iPer, iMax, nProc, skeletonFile);
	}
	else if(execType == "waitSlurm")
	{
		int numberJobs;
		stringstream cGen;
		cGen << argv[2];
		cGen >> numberJobs;
		waitSlurmFinish(numberJobs);

	}
	else if(execType == "buildCsvFile")
	{
		int size;
		string skeletonName;
		stringstream cGen;
		cGen << argv[2] << "  " << argv[3];
		cGen >> size >> skeletonName;
		buildCsvFile(size, skeletonName);
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


void waitSlurmFinish(int usedProc)
{
        while(true)
        {
#ifdef UNIX
                sleep(0.1);
#endif
                system("squeue > quee.txt");
                ifstream in_("quee.txt");
                string line;
                for(int i = 0; i < usedProc; i++)
                {
                        getline(in_,line);
                }
                in_.close();
                remove("quee.txt");

                if(line == "")
                        break;
        }
}


void buildCsvFile(int size, string skeletonName)
{
	system("cat *.csv > allIsomers.csv");
	string csvName = "allIsomers";

	ifstream csvFile1_((csvName + ".csv").c_str());
	ifstream skeleton_(skeletonName.c_str());
	ofstream csvFileNumber_((csvName + "-number.csv").c_str());

	string line1;
	int k1 = 1;
	while(getline(csvFile1_,line1))
	{
		if(line1 == "")
		{
			csvFileNumber_ << endl;
			continue;
		}


		string aux1, aux2;
		vector<int> permutCsv(size);
		stringstream lineCsv;
		lineCsv << line1;
		lineCsv >> aux1 >> aux2;
		for(int i = 0; i < size; i++)
			lineCsv >> permutCsv[i];

		string line2;
		while(getline(skeleton_,line2))
		{
			if(line2 == "")
				continue;

			vector<int> permutSkeleton(size);
			stringstream lineSkel;
			lineSkel << line2;
			for(int i = 0; i < size; i++)
				lineSkel >> permutSkeleton[i];
			
			if(permutSkeleton == permutCsv)
			{
				csvFileNumber_ << aux1 << "  " 
				<< aux2 << " " << k1 << " ; ";
				for(int i = 0; i < size; i++)
					csvFileNumber_ << permutCsv[i] << "  ";

				csvFileNumber_ << endl;
				k1++;
				break;
			}
			k1++;
		}
	}
	csvFileNumber_.close();
	skeleton_.close();
	csvFile1_.close();

	system("cat *weights* > allWeights.txt");
	vector<int> allWeights;
	ifstream weightsFile_("allWeights.txt");
	if(weightsFile_.is_open())
	{
		string line;
		while(getline(weightsFile_,line))
		{
			stringstream readLine;
			readLine << line;
			int number;
			readLine >> number;
			allWeights.push_back(number);
		}
		sort(allWeights.begin(),allWeights.end());
	}

	ifstream csvNumber_((csvName + "-number.csv").c_str());
	ofstream csvFinal_((csvName + "-final.csv").c_str());
	int kAll = 0;
	if(allWeights.size() == 0)
		allWeights.push_back(-1);

	string line;
	while(getline(csvNumber_,line))
	{
		if(line == "")
		{
			csvFinal_ << endl;
			continue;
		}

		string dum1, dum2;
		vector<int> permutCsv(size);
		int currentWeight, permutIndex;
		stringstream csvLine;
		csvLine << line;
		csvLine >> currentWeight >> dum1 >> permutIndex >> dum2;
		for(int i = 0; i < size; i++)
			csvLine >> permutCsv[i];

		while(permutIndex == allWeights[kAll])
		{
			currentWeight++;
			kAll++;
		}

		csvFinal_ << currentWeight << "  ;  " << permutIndex << "  ;  ";
		for(int i = 0; i < size; i++)
			csvFinal_ << permutCsv[i] << "  ";
		csvFinal_ << endl;

	}

}


void changeNameOfFiles(string responseName)
{
	ifstream response_(responseName.c_str());
	string line;
	ofstream counting_("counting.csv");
	while (!response_.eof())
	{
		getline(response_, line);
		string combination;
		stringstream convert;
		convert << line;
		convert >> combination;
		AllMolecularFormulas allMol_;
		vector< vector<int> > combinationCode = allMol_.stringToNumber(combination);
		string newCombinationName = allMol_.newCodeToString(combinationCode);
		int systemSize = 0;
		for (size_t i = 0; i < combinationCode.size(); i++)
		{
			for (size_t j = 0; j < combinationCode[i].size(); j++)
			{
				if(i > 0)
					systemSize += 2 * combinationCode[i][j];
				else
					systemSize += combinationCode[i][j];
			}
		}
		string geomName = sizeToGeometryCode(systemSize);
		ofstream newFile_((geomName + "-" + newCombinationName + ".csv").c_str());
		ifstream typesFile_((combination + "---atomTypes.txt").c_str());
		string typeLine;
		getline(typesFile_, typeLine);
		newFile_ << typeLine << endl;
		ifstream isomerFile_((("final-") + combination).c_str());
		string isomerLine;

		// sum up weights
		int nWeights = 0;
		while (!isomerFile_.eof())
		{
			getline(isomerFile_, isomerLine);
			if (isomerLine == "")
				continue;
			stringstream convertWeights;
			int auxWeight;
			convertWeights << isomerLine;
			convertWeights >> auxWeight;
			nWeights += (auxWeight + 1);
		}
		isomerFile_.close();
		isomerFile_.open((("final-") + combination).c_str());

		int totalChiral = 0;
		int totalAchiral = 0;

		bool chiral;
		CauchyIndex cauchy_(systemSize);
		while (!isomerFile_.eof())
		{
			getline(isomerFile_, isomerLine);
			if (isomerLine == "")
			{
				newFile_ << endl;
				continue;
			}

			//check next to see if chiral
			string isomerLine2;
			if (isomerFile_.eof())
				chiral = false;
			else
			{
				getline(isomerFile_, isomerLine2);
				if (isomerLine2 == "")
					chiral = false;
				else
					chiral = true;
			}
			stringstream convertLine1;
			int auxWeight;
			string dummy1, dummy2, dummy3;
			convertLine1 << isomerLine;
			convertLine1 >> auxWeight;
			convertLine1 >> dummy1 >> dummy2 >> dummy3;
			vector<int> notation1(systemSize);
			for (int i = 0; i < systemSize; i++)
				convertLine1 >> notation1[i];
			
			newFile_ << (auxWeight + 1) << " ; "
				<< "{" << geomName << " "
				<< "[" << notation1[0];
			for (size_t i = 1; i < notation1.size(); i++)
				newFile_ << " " << notation1[i];
			newFile_ << "] ";

			newFile_ << fixed << setprecision(2)
				<< (100.0e0 * (double)(auxWeight + 1) / (double)nWeights) << "% ";
			if (!chiral)
			{
				newFile_ << "A}" << endl << endl;
				totalAchiral++;
			}
			else
			{
				totalChiral++;
				newFile_ << "C}" << endl;
				stringstream convertLine2;
				int auxWeight2;
				convertLine2 << isomerLine2;
				convertLine2 >> auxWeight2;
				convertLine2 >> dummy1 >> dummy2 >> dummy3;
				vector<int> notation2(systemSize);
				for (int i = 0; i < systemSize; i++)
					convertLine2 >> notation2[i];

				newFile_ << (auxWeight2 + 1) << " ; "
					<< "{" << geomName << " "
					<< "[" << notation2[0];
				for (size_t i = 1; i < notation2.size(); i++)
					newFile_ << " " << notation2[i];
				newFile_ << "] ";
				newFile_ << fixed << setprecision(2)
					<< (100.0e0 * (double)(auxWeight + 1) / (double)nWeights) << "% ";
				newFile_ << "C}" << endl;
			}

		}
		counting_ << newCombinationName << " ; " << totalChiral << "; " << totalAchiral << endl;
	}
}


string sizeToGeometryCode(int size)
{
	switch (size)
	{
	case 4:
		return "T-4";
		break;
	case 5:
		return "TBPY-5";
		break;
	case 6:
		return "OC-6";
		break;
	case 7:
		return "COC-7";
		break;
	case 8:
		return "SAPR-8";
		break;
	case 9:
		return "JTCTPR-9";
		break;
	case 10:
		return "JMBIC-10";
		break;
	default:
		cout << "size not found" << endl;
		exit(1);
	}
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
