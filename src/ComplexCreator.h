#ifndef COMPLEXCREATOR_H
#define COMPLEXCREATOR_H

#include <vector>
#include <string>

#include "Ligand.h"

class ComplexCreator
{
public:
	ComplexCreator(std::vector<Ligand> &allLigands_in);

	ComplexCreator(
		std::vector<Ligand> &allLigands,
		int maxChelation_in,
		int saMaxIterations_in,
		double maxAlfaAngle_in,
		double maxBetaAngle_in,
		double saTemperatureUpdate_in,
		double saInitialTemperature_in,
		double saAcceptance
		);

	~ComplexCreator();

	void calculateAllAngles(int nPoints);

	bool start();

	bool start(std::vector<int> & ligandsPermutation);

	std::vector<CoordXYZ> simulatedAnnealing();

	int getIterationOfLowestFit() { return finalI; }

	std::vector<Ligand> getLigandsCreated() const;

private:
	int maxChelation;
	double maxAlfaAngle;
	double maxBetaAngle;
	double saTemperatureUpdate;
	double saInitialTemperature;
	double saAcceptance;
	int saMaxIterations;
	int finalI;

	std::vector<Ligand> & allLigands;

	//start
	int orderAllLigands(); //bigger teeth first

	std::vector<double> getPoints(int totalChelation);

	void setInitialPosition(const std::vector<double> & points, std::vector<int> & ligandsPermutation);

	void setInitialPositionCauchy(const std::vector<double> & points, std::vector<int> & ligandsPermutation);

	std::vector<int> findGoodPoint(
		int chelation, 
		const std::vector<double> & points,
		std::vector<bool> & pointsTaken);

	int closestPoint(
		double x, double y, double z,
		const std::vector<double> &points,
		std::vector<bool> &pointsTaken
		);

	//optimizing
	double calculateAllFit(std::vector<Ligand> & ligands);

	void perturbOperations(std::vector<Ligand> & ligands);

	void printAllAtoms(std::string xyzName, std::vector<Ligand> & ligands);
};


#endif


