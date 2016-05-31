#ifndef ROOTMEANSQUAREDEVIATION_H
#define ROOTMEANSQUAREDEVIATION_H

#include <string>
#include <vector>

#include "Coordstructs.h"

class RootMeanSquareDeviation
{
public:
	RootMeanSquareDeviation();
	
	~RootMeanSquareDeviation();

	double rmsd(std::string file1, std::string file2);

	double rmsOverlay(std::string molName1, std::string molName2);

	double rmsOverlay(std::vector<CoordXYZ> & mol1, std::vector<CoordXYZ> & mol2);

	std::vector<CoordXYZ> readCoord(std::string fName);

private:
	std::vector<double> readPoint(std::string fName, int format = 1);
	std::vector<double> rotateToZ0(const std::vector<double> &point);
	std::vector<double> rotateToPlane(const std::vector<double> &point);
	std::vector<double> mirrorY(const std::vector<double> &point);
	void printXyz(std::string fName, const std::vector<double> &points);
	void printXyz(std::string fName, std::vector<CoordXYZ> &mol);
	void printXyzSuperpositions(std::string fName, std::vector<CoordXYZ> &mol1, std::vector<CoordXYZ> &mol2);
	double rms(std::vector<CoordXYZ> &mol1, std::vector<CoordXYZ> &mol2);
	void rotateMol(std::vector<CoordXYZ> &mol, double x, double y, double z, double angle);
};


#endif