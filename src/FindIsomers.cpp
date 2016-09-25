#include "FindIsomers.h"

#include <vector>
#include <string>

#include "BuildComplex.h"
#include "Coordstructs.h"

using namespace std;

FindIsomers::~FindIsomers(){}

FindIsomers::FindIsomers()
{
	// B e C are different atoms
	// limpar os arquivos dos ligantes e fazer as permutacoes necessarias

	BuildComplex bc_;
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA();
	int permutationsNumber = bc_.getLigandsPermutation().size();

}
