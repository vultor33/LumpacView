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
	log_ << rmsd_.rmsd("hardin4.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin5.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin6.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin7.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin8.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin9.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin10.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin11.txt", "repulsion6.txt") << endl;
	log_ << rmsd_.rmsd("hardin12.txt", "repulsion6.txt") << endl;
	log_.close();
}

PointAnalysis::~PointAnalysis(){}

