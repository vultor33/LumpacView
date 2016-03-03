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
		int maxChelation_in,
		double stretchDistance_in,
		double maxAlfaAngle_in,
		double maxBetaAngle_in,
		int saMaxIterations_in,
		double saTemperatureUpdate_in
		);
	~ComplexCreator();

	bool start();

	std::vector<CoordXYZ> simulatedAnnealing();

private:
	int maxChelation;
	double stretchDistance;
	double maxAlfaAngle;
	double maxBetaAngle;
	double saTemperatureUpdate;
	int saMaxIterations;

	std::vector<Ligand> & allLigands;

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
