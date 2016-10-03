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

	std::vector<CoordXYZ> build(std::string ligandName, int coordination, int charge, std::vector<std::string> options, std::string mopacExecPath, int & nLigandAtoms);

	std::vector<CoordXYZ> build(ReadInput & readInp_);

	void makeComplexOptimizingInMopac(std::string ligandName, int coordination, int charge, std::vector<std::string> options, std::string mopacExecPath);

	std::vector<Ligand> BuildComplex::assembleComplexWithoutSA(std::vector<int> & ligandsPermutation = std::vector<int>(), std::vector< std::string > & inputInformations = std::vector< std::string >());

	void fitSA();

	std::vector<int> getLigandsPermutation() { return actualLigandPermutation; }

	int getLigandsNumber();

private:
	ReadInput BuildComplex::activateReadInput(std::vector< std::string > & inputInformations = std::vector< std::string >());

	std::vector<int> actualLigandPermutation;

	bool ReadLumpacViewInput(ReadInput & readInp_);

	bool buildLigands(ReadInput & readInp_); // set ligands directions X1 and X2.

	bool constructComplex(ReadInput & readInp_, const AdjustSaParameters & saParameters_, std::vector<CoordXYZ> &allAtoms); // place ligands on the sphere and optimize with simulated annealing.
	
	void runMopac(ReadInput & readInp_, std::vector<CoordXYZ> & allAtoms); // run optimization and frequency with previous coordinates.

	void setClandH2oToCompleteCoordination(Ligand & Cl_, Ligand & H2o_);

	bool optimize(
		std::string mopacExePath,
		std::vector<CoordXYZ> & allAtoms,
		std::vector<std::string> & options);

	bool optimize(
		std::string mopacExePath,
		std::vector<CoordXYZ> & allAtoms,
		std::vector<std::string> & options,
		std::vector<MopacParams> & params);

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