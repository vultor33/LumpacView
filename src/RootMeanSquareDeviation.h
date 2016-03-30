#ifndef ROOTMEANSQUAREDEVIATION_H
#define ROOTMEANSQUAREDEVIATION_H

#include <string>
#include <vector>

class RootMeanSquareDeviation
{
public:
	RootMeanSquareDeviation();
	
	~RootMeanSquareDeviation();

	double rmsd(std::string file1, std::string file2);

private:
	std::vector<double> readPoint(std::string fName);
	std::vector<double> rotateToZ0(const std::vector<double> &point);
	std::vector<double> rotateToPlane(const std::vector<double> &point);
	std::vector<double> mirrorY(const std::vector<double> &point);
	void printXyz(std::string fName, const std::vector<double> &points);




};


#endif