#ifndef GEOMETRIES_H
#define GEOMETRIES_H

#include <vector>
#include <string>

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
	
	void selectGeometrySymmetries(
		int select,
		std::vector< std::vector<int> > &allReflections);

	std::string selectGeometrySymmetriesFlag(
		int select,
		int iSymmetry,
		int symmetryType);

	std::string sizeToGeometryCode(int size);



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

	std::vector<double> geometry9TCTPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9CSAPR(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9MFF(
        	std::vector<CoordXYZ> &mol0,
        	double & cutAngle,
        	std::vector<int> &reflectionOperation);

	std::vector<double> geometry10PointSphere(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry10TD(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry10JSPC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry10JBCSAPR(
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




	///////////////////////////////////
	////////// REFLECTIONS ////////////
	///////////////////////////////////
	void geometry5SPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry5SPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry5TBPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry5TBPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);


	void geometry6OCotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry6OCSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry6TPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry6TPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);


	void geometry7COCotherSymmetries(std::vector< std::vector<int> > &allReflections);
	std::string geometry7COCSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7CTPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7CTPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8BTPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8BTPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);


	void geometry9CSAPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9CSAPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9MFFotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9MFFSymmetryFlags(
		int iSymmetry,
		int symmetryType);


	AuxMath auxMath_;


};


#endif
