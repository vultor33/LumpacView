#ifndef BUILDCOMPLEX_H
#define BUILDCOMPLEX_H

#include <vector>
#include <string>

#include "AdjustSaParameters.h"
#include "ReadInput.h"
#include "Coordstructs.h"
#include "Ligand.h"

class BuildComplex
{
public:
	BuildComplex();

	~BuildComplex();

	std::vector<CoordXYZ> build();

	std::vector<CoordXYZ> build(ReadInput & readInp_);

	std::vector<CoordXYZ> buildCompletingWithWater(std::string ligandName, int coordination, int charge, std::vector<std::string> options, std::string mopacExecPath, int & nLigandAtoms);

	void makeComplexOptimizingInMopac(std::string ligandName, int coordination, int charge, std::vector<std::string> options, std::string mopacExecPath);
	
	std::vector<Ligand> assembleComplexWithoutSA();

	std::vector<Ligand> assembleComplexWithoutSA(
		std::vector<int> & ligandsPermutation);

	std::vector<Ligand> assembleComplexWithoutSA(
		std::vector<int> & ligandsPermutation,
		std::vector< std::string > & inputInformations);

	std::vector<Ligand> assembleComplexWithoutSACauchy(
		std::vector<int> & ligandsPermutation, 
		std::vector< std::string > & inputInformations,
		std::string flagsFile);
		
	void runMopacAndPrint(std::vector< std::string > options, std::string mopacExecPath, std::vector<CoordXYZ> & allAtoms); // run optimization and frequency with previous coordinates.

	void fitSA();

	std::vector<int> getLigandsPermutation() { return actualLigandPermutation; }

	int getLigandsNumber();

private:
	ReadInput activateReadInput();

	ReadInput activateReadInput(std::vector< std::string > & inputInformations);

	std::vector<int> actualLigandPermutation;

	bool ReadLumpacViewInput(ReadInput & readInp_);

	bool buildLigands(ReadInput & readInp_); // set ligands directions X1 and X2.

	bool constructComplex(ReadInput & readInp_, const AdjustSaParameters & saParameters_, std::vector<CoordXYZ> &allAtoms); // place ligands on the sphere and optimize with simulated annealing.
	
	void runMopac(ReadInput & readInp_, std::vector<CoordXYZ> & allAtoms); // run optimization and frequency with previous coordinates.

	void runMopac(std::vector< std::string > options, std::string mopacExecPath, std::vector<CoordXYZ> & allAtoms); // run optimization and frequency with previous coordinates.

	void runMopacAndGetCoordinates(std::vector< std::string > options, std::string mopacExecPath, std::vector<CoordXYZ> & allAtoms); 

	void setClandH2oToCompleteCoordination(Ligand & Cl_, Ligand & H2o_);

	bool optimize(
		std::string mopacExePath,
		std::vector<CoordXYZ> & allAtoms,
		std::vector<std::string> & options,
		std::vector<MopacParams> params = std::vector<MopacParams>());

	void createParamsFile(
		std::string paramsName, 
		std::vector<MopacParams>& params);
};

#endif


/*
EXEMPLOS DE USOS
=======================================================================================================================================
std::vector<Ligand> BuildComplex::assembleComplexWithoutSA(std::vector<int> & ligandsPermutation = std::vector<int>(), std::vector< std::string > & inputInformations = std::vector< std::string >());
->
vector< string > inputInformations(size);
inputInformations[0] = "metalName";
inputInformations[1] = "metalParams";
inputInformations[2] = "agua";
inputInformations[3] = "agua";
inputInformations[4] = "agua";
=======================================================================================================================================
*/


/* MAKING COMPLEX

-> FICAR ATENTO A CAPTACAO DO NOVO LIGANTE

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
