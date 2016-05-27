#ifndef COORDSTRUCTS_H
#define COORDSTRUCTS_H

#include <string>
#include <vector>

struct CoordXYZ
{
	std::string atomlabel;
	double x, y, z;
};

struct MopacParams
{
	std::string paramName;
	double paramValue;
};

/*
struct Ligand
{
	std::vector< CoordXYZ > coord;
	int x1;
	int x2;

	// r, teta, fi
	// alfa de rotacao solida

	double teta;
	double ligx, ligy, ligz;

};
*/


#endif