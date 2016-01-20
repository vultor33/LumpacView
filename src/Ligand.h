#ifndef LIGAND_H
#define LIGAND_H

#include "Coordstructs.h"

class Ligand
{
public:
	Ligand();
	~Ligand();

	void initializeLigand(std::vector<CoordXYZ> &coord_in);

private:
	std::vector<CoordXYZ> coord;
	int x1; // donor atom
	int x2; // ligand orientation vactor

};


#endif
