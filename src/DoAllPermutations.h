#ifndef DOALLPERMUTATIONS_H
#define DOALLPERMUTATIONS_H

#include <string>
#include <vector>

class DoAllPermutations
{
public:
	DoAllPermutations();

	~DoAllPermutations();

	void calculateAll(
		std::string methodOptimize,
		std::string methodCosmo,
		std::string projectName,
		std::string filePermutations,
		int nPermutations);

	void analysis(std::string crystalFile, std::string isomerFile);

private:
	std::vector< std::vector<int> > permutations;

	void getAllPermutations(std::string filePermutations, int nPermutations);

	void buildComplexWithAPermutation(
		std::vector<int> permutation,
		std::string methodOptimize,
		std::string methodCosmo,
		std::string projectName
		);

	std::string permutationToString(std::vector<int>& permutation);

};

/* EXEMPLO DE USO
DoAllPermutations doall_;
string methodOptimize = " RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
string methodCosmo = " EPS=78.4 1SCF PRECISE NOINTER XYZ T=10D + \n NOLOG GEO-OK SCFCRT=1.D-10";
int nPermutations = 7;
string projectName = "jalniu-RM1";
string filePermutations = "BBm1m1m1.xyz";
doall_.calculateAll(
	methodOptimize,
	methodCosmo,
	projectName,
	filePermutations,
	nPermutations);
*/


#endif