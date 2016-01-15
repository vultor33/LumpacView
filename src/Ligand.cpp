#include "Ligand.h"

#include "Coordstructs.h"

using namespace std;

Ligand::Ligand() {}

Ligand::~Ligand() {}

void Ligand::initializeLigand(vector<CoordXYZ> &coord_in, int x1_in, int x2_in)
{
	coord = coord_in;
	x1 = x1_in;
	x2 = x2_in;
}

