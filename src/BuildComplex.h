#ifndef BUILDCOMPLEX_H
#define BUILDCOMPLEX_H

#include <vector>

#include "AdjustSaParameters.h"
#include "ReadInput.h"
#include "Coordstructs.h"

class BuildComplex
{
public:
	BuildComplex();
	~BuildComplex();

	void build();

	void build(std::string ligandName, int coordination, int charge, std::vector<std::string> options);

	void build(ReadInput & readInp_);

	void checkIfIsSameIsomer(std::string xRayName);

	void fitSA();


private:
	bool ReadLumpacViewInput(ReadInput & readInp_);
	bool buildLigands(ReadInput & readInp_); // set ligands directions X1 and X2.
	bool constructComplex(ReadInput & readInp_, const AdjustSaParameters & saParameters_, std::vector<CoordXYZ> &allAtoms); // place ligands on the sphere and optimize with simulated annealing.
	void runMopac(ReadInput & readInp_, std::vector<CoordXYZ> & allAtoms); // run optimization and frequency with previous coordinates.

	void setClandH2oToCompleteCoordination(Ligand & Cl_, Ligand & H2o_);


};

#endif