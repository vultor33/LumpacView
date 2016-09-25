#include "FindIsomers.h"

#include <vector>
#include <string>

#include "BuildComplex.h"
#include "Coordstructs.h"

using namespace std;

FindIsomers::~FindIsomers(){}

FindIsomers::FindIsomers()
{
	BuildComplex bc_;
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA();
	int permutationsNumber = bc_.getLigandsPermutation().size();

}
