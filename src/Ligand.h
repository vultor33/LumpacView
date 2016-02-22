#ifndef LIGAND_H
#define LIGAND_H

#include <vector>
#include <string>

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

	bool initializeLigand2();

private:
	//data
	std::vector<CoordXYZ> coord;
	std::string titleInfo;	
	int chelation;
	CoordXYZ X1; //coordination center
	CoordXYZ X2; //vector - leaving from X1

	bool getInfoFromTitle();
	void printXyzLigandDirection(); //debug purpose

	bool calculateMonodentate();
	bool calculateBidentate();
	bool calculateTridentate();


};

#endif



