#ifndef COMPLEXCREATOR_H
#define COMPLEXCREATOR_H

#include <vector>
#include <string>

#include "Ligand.h"

class ComplexCreator
{
public:
	ComplexCreator(
		std::vector<Ligand> &allLigands,
		std::string metalName_in,
		std::string metalParams_in,
		std::string projectName_in
		);
	~ComplexCreator();

	void start();


private:
	std::vector<Ligand> & allLigands;
	std::string metalName;
	std::string metalParams;
	std::string projectName;


};



#endif
