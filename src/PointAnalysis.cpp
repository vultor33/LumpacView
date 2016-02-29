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
	log_ << rmsd_.rmsd("repulsion5.txt", "fred5.txt") << endl;
	log_ << rmsd_.rmsd("repulsion6.txt", "fred6.txt") << endl;
	log_ << rmsd_.rmsd("repulsion7.txt", "fred7.txt") << endl;
	log_ << rmsd_.rmsd("repulsion8.txt", "fred8.txt") << endl;
	log_ << rmsd_.rmsd("repulsion9.txt", "fred9.txt") << endl;
	log_ << rmsd_.rmsd("repulsion10.txt", "fred10.txt") << endl;
	log_ << rmsd_.rmsd("repulsion11.txt", "fred11.txt") << endl;
	log_ << rmsd_.rmsd("repulsion12.txt", "fred12.txt") << endl;
	log_.close();
}

PointAnalysis::~PointAnalysis(){}

