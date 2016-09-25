#include "FindIsomers.h"

#include <vector>
#include <string>

#include "BuildComplex.h"
#include "Coordstructs.h"

using namespace std;

FindIsomers::~FindIsomers(){}

/*

MUDANCA IMPORTANTE -> NO BIDENTADO DE 2 ATOMOS QUALQUER LADO E PERPENDICULAR
                      NO TRIDENTADO DE 3 ATOMOS BASTA CALCULAR O PRODUTO VETORIAL E PRONTO

*/

FindIsomers::FindIsomers()
{
	// B e C are different atoms

	BuildComplex bc_;
	vector<Ligand> allAtomsOriginal = bc_.assembleComplexWithoutSA();
	int permutationsNumber = bc_.getLigandsPermutation().size();

}
