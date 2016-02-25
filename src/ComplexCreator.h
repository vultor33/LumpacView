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

	void simulatedAnnealing();

private:
	const int maxChelation = 10;
	const double stretchDistance = 3;
	const double maxAlfaAngle = 3.14159265358979323846 / 90.0e0;
	const double maxBetaAngle = 3.14159265358979323846 / 90.0e0;
	const int saMaxIterations = 1000;
	
	//data
	std::vector<Ligand> & allLigands;
	std::string metalName;
	std::string metalParams;
	std::string projectName;

	//start
	int orderAllLigands(); //bigger teeth first
	std::vector<double> getPoints(int totalChelation);
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

	//optimizing
	double calculateAllFit(std::vector<Ligand> & ligands);

	void perturbOperations(std::vector<Ligand> & ligands);

	void printAllAtoms(std::vector<Ligand> & ligands);
};


#endif
