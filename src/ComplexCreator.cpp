#include "ComplexCreator.h"

#include <vector>

#include "Ligand.h"

using namespace std;

ComplexCreator::ComplexCreator(
	vector<Ligand> &allLigands_in,
	string metalName_in,
	string metalParams_in,
	string projectName_in
	)
	:allLigands(allLigands_in)
{
	metalName = metalName_in;
	metalParams = metalName_in;
	projectName = projectName_in;
}

ComplexCreator::~ComplexCreator(){}

void ComplexCreator::start()
{
	size_t size = allLigands.size();
	int sumChelation = 0;
	for (size_t i = 0; i < size; i++)
	{
		sumChelation += allLigands[i].getChelation();
	}



}

// posso seguir na ideia do Simas.
// posso estabelecer a situacao pela construcao dos pontos.
















