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
		int saMaxIterations_in,
		double maxAlfaAngle_in,
		double maxBetaAngle_in,
		double saTemperatureUpdate_in,
		double saInitialTemperature_in,
		double saAcceptance
		);
	~ComplexCreator();

	bool start();

	std::vector<CoordXYZ> simulatedAnnealing();

#ifdef _FITSA
	double finalFit = 466;
	int finalI = 5000;
#endif

private:
	int maxChelation;
	double stretchDistance;
	double maxAlfaAngle;
	double maxBetaAngle;
	double saTemperatureUpdate;
	double saInitialTemperature;
	double saAcceptance;
	int saMaxIterations;

	std::vector<Ligand> & allLigands;

	//start
	int orderAllLigands(); //bigger teeth first

	std::vector<double> getPoints(int totalChelation);

	void setInitialPosition(const std::vector<double> & points);

	std::vector<int> findGoodPoint(
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
