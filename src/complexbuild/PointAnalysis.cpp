#include "PointAnalysis.h"

#include <fstream>

#include "RootMeanSquareDeviation.h"
#include "AuxMath.h"

using namespace std;

PointAnalysis::PointAnalysis()
{
	RootMeanSquareDeviation rmsd_;
	
	ofstream log_("log.txt");

	log_ << "RMSD" << endl;
//	log_ << rmsd_.oldRmsdToAddressPoints("hardin4.txt", "repulsion6.txt") << endl;
	log_.close();
}

PointAnalysis::~PointAnalysis(){}

