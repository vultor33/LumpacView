#ifndef LIGAND_H
#define LIGAND_H

#include "Coordstructs.h"

class Ligand
{
public:
	Ligand();
	~Ligand();

	void setLigandCoordinates(std::vector<CoordXYZ> &coord_in);
	bool initializeLigand();

private:
	std::vector<CoordXYZ> coord;
	int chelation; // mono, bi or tri

	void printXyzLigandDirection(); //debug purpose

};

#endif



