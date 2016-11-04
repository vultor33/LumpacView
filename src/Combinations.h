#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include <vector>
#include <string>

class Combinations
{
public:
	Combinations();

	~Combinations();

	void doAllCombinations();

	std::vector< std::vector<int> > stringToNumber(std::string entryString);

private:

	int codeToType(std::string code);

	void addEqual(
		int codeNumber,
		std::vector<int> & typeCode1,
		std::vector<int> & typeCode2,
		std::vector<int> & typeCode3
		);

	void addDifferent(
		int codeNumber,
		std::vector<int> & typeCode1,
		std::vector<int> & typeCode2,
		std::vector<int> & typeCode3
		);






};


#endif
