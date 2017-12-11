#ifndef GEOMETRIES_H
#define GEOMETRIES_H

#include <vector>

#include "Coordstructs.h"
#include "AuxMath.h"

class Geometries
{
public:
	Geometries();
	~Geometries();
	
	std::vector<double> selectGeometry(
		int select,
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);


private:
	std::vector<double> geometry4Tetrahedron(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry4Square(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry5TBPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry5SPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry5VOC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry6OC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry6TPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry7COC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry7PBPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry7CTPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8SAPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8TDD(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8BTPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8HBPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8CU(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9JTCTPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9CSAPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry10JMBIC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry11JCPAPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry12IC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);




	AuxMath auxMath_;


};


#endif