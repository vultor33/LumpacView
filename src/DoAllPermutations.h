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
		std::string projectName);

private:
	std::vector< std::vector<int> > permutations;

	void setAllPermutationsManually(int nComplex);

	void buildComplexWithAPermutation(
		std::vector<int> permutation,
		std::string methodOptimize,
		std::string methodCosmo,
		std::string projectName
		);



};


#endif