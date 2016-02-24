#ifndef LIGAND_H
#define LIGAND_H

#include <vector>
#include <string>
#include <fstream>

#include "Coordstructs.h"

class Ligand
{
public:
	Ligand();
	~Ligand();

	void setLigandCoordinates(
		std::vector<CoordXYZ> &coord_in, 
		std::string titleInfo_in);
	
	bool initializeLigand();

	inline int getChelation() { return chelation; }

	void translateLigand(double x, double y, double z);

	void rotateToCenter();

	void rotateOnX1(double vx, double vy, double vz, double ang);

	void printLigand(std::ofstream &out);

	inline int getNatoms() { return (int)coord.size(); }

private:
	//data
	std::vector<CoordXYZ> coord;
	std::string titleInfo;	
	int chelation;
	CoordXYZ X1; //coordination center
	CoordXYZ X2; //vector - leaving from X1

	bool getInfoFromTitle();
	bool calculateMonodentate();
	bool calculateBidentate();
	bool calculateTridentate();

	void printXyzLigandDirection(std::string inputName = "xyzTeste.xyz"); //debug purpose

};

#endif



