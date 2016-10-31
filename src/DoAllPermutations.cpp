#include "DoAllPermutations.h"

#include <vector>
#include <string>
#include <sstream>

#include "Coordstructs.h"
#include "FindIsomers.h"
#include "BuildComplex.h"

using namespace std;

DoAllPermutations::DoAllPermutations() {}

DoAllPermutations::~DoAllPermutations(){}

void DoAllPermutations::calculateAll(string methodOptimize, string methodCosmo, string projectName)
{



}

void DoAllPermutations::buildComplexWithAPermutation(
	std::vector<int> permutation,
	std::string methodOptimize,
	std::string methodCosmo,
	std::string projectName
	)
{
	FindIsomers fd_;
	vector<string> options;
	string mopacExecPath;
	vector<CoordXYZ> allAtoms = fd_.buildComplexWithSelectedIsomer(
		permutation,
		projectName,
		methodOptimize,
		options,
		mopacExecPath);
	BuildComplex bc_;
	options[1] = projectName + "-water";
	options[2] = methodCosmo;
	options[3] = "";
	options[4] = "";
	bc_.runMopacAndPrint(options, mopacExecPath, allAtoms);
}

void DoAllPermutations::setAllPermutationsManually(int nComplex)
{

	if (nComplex == 0)
	{
		permutations.resize(12);
		for (size_t i = 0; i < 12; i++)
			permutations[i].resize(7);

		permutations[0][0] = 0;


	}
	else
		exit(1);
}





