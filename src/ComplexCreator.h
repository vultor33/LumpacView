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
	const double stretchDistance = 5;

	std::vector<Ligand> & allLigands;
	std::string metalName;
	std::string metalParams;
	std::string projectName;

	int orderAllLigands(); //bigger teeth first
	std::vector<double> getPoints(int totalChelation);
	std::vector<double> arrayToVector(const double * array_in, size_t size);

	void setInitialPosition(const std::vector<double> & points);
	
	std::vector<double> findGoodPoint(
		int chelation, 
		const std::vector<double> & points,
		std::vector<bool> & pointsTaken);

	int closestPoint(
		double x, double y, double z,
		const std::vector<double> &points,
		std::vector<bool> &pointsTaken
		);

	void stretchPoints(std::vector<double> &points);

	void printAllAtoms();


};


#endif
