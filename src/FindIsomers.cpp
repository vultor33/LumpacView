#include "FindIsomers.h"

#include <vector>
#include <string>

#include "BuildComplex.h"
#include "Coordstructs.h"

using namespace std;

FindIsomers::~FindIsomers(){}

FindIsomers::FindIsomers()
{
	// H, B and C are different atoms

	BuildComplex bc_;
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA();
	int permutationsNumber = bc_.getLigandsPermutation().size();

}
