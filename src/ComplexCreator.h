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

	bool start();


private:
	const int maxChelation = 10;

	std::vector<Ligand> & allLigands;
	std::string metalName;
	std::string metalParams;
	std::string projectName;

	int orderAllLigands(); //bigger teeth first
	std::vector<double> getPoints(int totalChelation);
	std::vector<double> arrayToVector(const double * array_in, size_t size);

};



#endif
