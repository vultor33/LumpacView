#include "Fitness.h"

#include <vector>

#include "AuxMath.h"
#include "Coordstructs.h"

using namespace std;

Fitness::Fitness()
{
	fitType = 0;
}

Fitness::~Fitness() {}

double Fitness::calculateFit(
	vector<CoordXYZ> & ligand1,
	vector<CoordXYZ> & ligand2)
{
	double fit = 0.0e0;
	double r;
	AuxMath auxMath_;
	if (fitType == 0)
	{
		//maximizing distance
		for (size_t i = 0; i < ligand1.size(); i++)
		{
			for (size_t j = 0; j < ligand2.size(); j++)
			{
				r = auxMath_.norm(
					ligand1[i].x - ligand2[j].x,
					ligand1[i].y - ligand2[j].y,
					ligand1[i].z - ligand2[j].z);
				fit += 1 / r;
			}
		}
	}

	return fit;
}

