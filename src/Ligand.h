#ifndef LIGAND_H
#define LIGAND_H

#include "Coordstructs.h"

class Ligand
{
public:
	Ligand();
	~Ligand();

	void initializeLigand(std::vector<CoordXYZ> &coord_in, int x1_in, int x2_in);

private:
	std::vector<CoordXYZ> coord;
	int x1; // centro da ligacao
	int x2; // orientacao

};


#endif
