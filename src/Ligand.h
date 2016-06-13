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

	void setLigandCoordinates(std::string ligandFileName);

	bool initializeLigand();

	inline int getChelation() { return chelation; }

	void translateLigand(double x, double y, double z);

	void rotateToCenter();

	void rotateOnX1(double vx, double vy, double vz, double ang);

	void rotateOnX2(double beta); //vector x1 x2

	void genericRotation(double vx, double vy, double vz, double ang);

	void placeLigandOnPoins(std::vector<int> &pLig,
		const std::vector<double> & points);

	void rotateOverReferencePoints(double angle = 0);

	double distanceX1ToPoint(double x, double y, double z);

	void printLigand(std::ofstream &out);

	void setNewCoordinates(std::vector<CoordXYZ> &newCoord);

	std::vector<CoordXYZ> getAllAtoms();

	inline int getNatoms() { return (int)coord.size(); }

	inline std::string getFormalName() { return formalName; }

private:
	//data
	std::vector<CoordXYZ> coord;
	std::vector<double> sphereReferencePoints; //xyz xyz xyz
	std::string titleInfo;
	double metalDistance;
	int chelation;
	std::string formalName;
	CoordXYZ X1; //coordination center
	CoordXYZ X2; //vector - leaving from X1

	bool getInfoFromTitle();
	bool calculateMonodentate();
	bool calculateBidentate();
	bool calculateTridentate();

	void printXyzLigandDirection(std::string inputName = "zxyzTeste.xyz"); //debug purpose
};

#endif



