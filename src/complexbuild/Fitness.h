#ifndef FITNESS_H
#define FITNESS_H

#include <vector>

#include "Coordstructs.h"

class Fitness
{
public:
	Fitness();
	~Fitness();

	double calculateFit(
		std::vector<CoordXYZ> & ligand1,
		std::vector<CoordXYZ> & ligand2);

private:
	int fitType;

};



#endif
