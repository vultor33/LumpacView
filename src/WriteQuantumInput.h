#ifndef WRITEQUANTUMINPUT_H
#define WRITEQUANTUMINPUT_H

#include <fstream>
#include <string>
#include <vector>

#include "Coordstructs.h"

class WriteQuantumInput
{
public:
	WriteQuantumInput(std::vector<std::string> options);
	~WriteQuantumInput();

	std::string createInput(
		std::vector<CoordXYZ> &coordinates,
		int indexAddedAtFinalOfInputToSaveItsName = 0);

	// GLOBAL PUBLIC FUNCTIONS
	void changeProjectName(std::string newProjectName) { projectName = newProjectName; }

	//MOPAC PUBLIC FUNCTIONS
	void changeMopacHeader(std::string newHeader) { mopacHeader = newHeader; }

private:
	//GLOBAL
	std::string type;
	void setInputProperties(std::vector<std::string> &options);
	std::string projectName;

	//MOPAC
	std::string mopacHeader;
	std::string externalParams;
	std::string centralMetal;
	void buildMopacInput(std::vector<CoordXYZ> &allAtoms, std::string inputName);

	//GAMESS
	std::vector<std::string> gamessHeader;
	std::vector<std::string> gamessAtomBasisFiles;
	std::vector<std::string> gamessEcpFiles;
	void readGamessAuxFiles();
	std::vector< std::vector<std::string> > atomBasis;
	std::vector< std::vector<std::string> > atomEcp;
	void buildGamessInput(std::vector<CoordXYZ> &coordinates, std::string inputName);
	int getChargeFromAtom(std::string atomName);

};

#endif

/*
mopac OPTIONS
0 -> mopac ou mopac2009
1 -> nome do projeto
2 -> cabecalho
3 -> parametros externos
4 -> metal central
*/