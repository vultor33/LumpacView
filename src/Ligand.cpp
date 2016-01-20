#include "Ligand.h"

#include "Coordstructs.h"

using namespace std;

Ligand::Ligand() {}

Ligand::~Ligand() {}

void Ligand::initializeLigand(vector<CoordXYZ> &coord_in)
{
	coord = coord_in;
}

