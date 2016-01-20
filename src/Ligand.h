#ifndef LIGAND_H
#define LIGAND_H

#include "Coordstructs.h"

class Ligand
{
public:
	Ligand();
	~Ligand();

	bool initializeLigand(std::vector<CoordXYZ> &coord_in);

private:
	std::vector<CoordXYZ> coord;
	int chelation; // mono, bi or tri

};

#endif



