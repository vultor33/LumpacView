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

	std::vector<double> geometry4SS(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry4vTBPY(
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

	std::vector<double> geometry5PP(
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

	std::vector<double> geometry6HP(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry6PPY(
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

	std::vector<double> geometry7HPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry7HP(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry7JETPY(
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

	std::vector<double> geometry8ETBPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8HPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8OP(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry8JGBF(
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

	std::vector<double> geometry9CCU(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9HH(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9OPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9EP(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9HBPY(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9JTC(
		std::vector<CoordXYZ> &mol0,
		double & cutAngle,
		std::vector<int> &reflectionOperation);

	std::vector<double> geometry9JTDIC(
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

	void geometry4TetrahedronotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry4TetrahedronSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry4SquareotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry4SquareSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry4SSotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry4SSSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry4vTBPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry4vTBPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

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

	void geometry5PPotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry5PPSymmetryFlags(
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

	void geometry6HPotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry6HPSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry6PPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry6PPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7COCotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7COCSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7CTPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7CTPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7PBPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7PBPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7HPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7HPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7HPotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7HPSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry7JETPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry7JETPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8BTPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8BTPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8HBPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8HBPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8SAPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8SAPRSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8TDDotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8TDDSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8CUotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8CUSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8ETBPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8ETBPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8HPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8HPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8OPotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8OPSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry8JGBFotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry8JGBFSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9TCTPRotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9TCTPRSymmetryFlags(
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

	void geometry9CCUotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9CCUSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9HHotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9HHSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9OPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9OPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9EPotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9EPSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9HBPYotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9HBPYSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9JTCotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9JTCSymmetryFlags(
		int iSymmetry,
		int symmetryType);

	void geometry9JTDICotherSymmetries(
		std::vector< std::vector<int> > &allReflections);
	std::string geometry9JTDICSymmetryFlags(
		int iSymmetry,
		int symmetryType);


	AuxMath auxMath_;


};


#endif
