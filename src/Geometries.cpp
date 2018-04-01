#include "Geometries.h"

#include <vector>
#include <iostream>
#include <stdlib.h>

#include "Coordstructs.h"
#include "AuxMath.h"

using namespace std;

Geometries::Geometries(){}

Geometries::~Geometries(){}

std::vector<double> Geometries::selectGeometry(
	int select,
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	switch (select)
	{
	case 40:
		return geometry4Tetrahedron(mol0, cutAngle, reflectionOperation);
		break;

	case 41:
		return geometry4Square(mol0, cutAngle, reflectionOperation);
		break;

	case 42:
		return geometry4SS(mol0, cutAngle, reflectionOperation);
		break;

	case 43:
		return geometry4vTBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 50:
		return geometry5TBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 51:
		return geometry5SPY(mol0, cutAngle, reflectionOperation);
		break;

	case 52:
		return geometry5VOC(mol0, cutAngle, reflectionOperation);
		break;

	case 53:
		return geometry5PP(mol0, cutAngle, reflectionOperation);
		break;

	case 60:
		return geometry6OC(mol0, cutAngle, reflectionOperation);
		break;

	case 61:
		return geometry6TPR(mol0, cutAngle, reflectionOperation);
		break;

	case 62:
		return geometry6HP(mol0, cutAngle, reflectionOperation);
		break;

	case 63:
		return geometry6PPY(mol0, cutAngle, reflectionOperation);
		break;

	case 70:
		return geometry7COC(mol0, cutAngle, reflectionOperation);
		break;

	case 71:
		return geometry7PBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 72:
		return geometry7CTPR(mol0, cutAngle, reflectionOperation);
		break;

	case 73:
		return geometry7HPY(mol0, cutAngle, reflectionOperation);
		break;

	case 74:
		return geometry7HP(mol0, cutAngle, reflectionOperation);
		break;

	case 75:
		return geometry7JETPY(mol0, cutAngle, reflectionOperation);
		break;

	case 80:
		return geometry8SAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 81:
		return geometry8TDD(mol0, cutAngle, reflectionOperation);
		break;

	case 82:
		return geometry8BTPR(mol0, cutAngle, reflectionOperation);
		break;

	case 83:
		return geometry8HBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 84:
		return geometry8CU(mol0, cutAngle, reflectionOperation);
		break;

	case 85:
		return geometry8ETBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 86:
		return geometry8HPY(mol0, cutAngle, reflectionOperation);
		break;

	case 87:
		return geometry8OP(mol0, cutAngle, reflectionOperation);
		break;

	case 88:
		return geometry8JGBF(mol0, cutAngle, reflectionOperation);
		break;

	case 90:
		return geometry9TCTPR(mol0, cutAngle, reflectionOperation);
		break;

	case 91:
		return geometry9CSAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 92:
		return geometry9MFF(mol0, cutAngle, reflectionOperation);
		break;

	case 93:
		return geometry9CCU(mol0, cutAngle, reflectionOperation);
		break;

	case 94:
		return geometry9HH(mol0, cutAngle, reflectionOperation);
		break;

	case 95:
		return geometry9OPY(mol0, cutAngle, reflectionOperation);
		break;

	case 96:
		return geometry9EP(mol0, cutAngle, reflectionOperation);
		break;

	case 97:
		return geometry9HBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 98:
		return geometry9JTC(mol0, cutAngle, reflectionOperation);
		break;

	case 99:
		return geometry9JTDIC(mol0, cutAngle, reflectionOperation);
		break;

	case 100:
		return geometry10PointSphere(mol0, cutAngle, reflectionOperation);
		break;

	case 101:
		return geometry10TD(mol0, cutAngle, reflectionOperation);
		break;

	case 102:
		return geometry10JSPC(mol0, cutAngle, reflectionOperation);
		break;

	case 103:
		return geometry10JBCSAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 110:
		return geometry11JCPAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 120:
		return geometry12IC(mol0, cutAngle, reflectionOperation);
		break;

	default:
		cout << "Geometry not found" << endl;
		exit(1);
		break;
	}
}

void Geometries::selectGeometrySymmetries(
	int select,
	std::vector< std::vector<int> > &allReflections)
{
	switch (select)
	{
	case 40:
		geometry4TetrahedronotherSymmetries(allReflections);
		break;

	case 41:
		geometry4SquareotherSymmetries(allReflections);
		break;

	case 42:
		geometry4SSotherSymmetries(allReflections);
		break;

	case 43:
		geometry4vTBPYotherSymmetries(allReflections);
		break;

	case 50:
		geometry5TBPYotherSymmetries(allReflections);
		break;

	case 51:
		geometry5SPYotherSymmetries(allReflections);
		break;

	case 52:
		//return geometry5VOC(mol0, cutAngle, reflectionOperation);
		break;

	case 53:
		geometry5PPotherSymmetries(allReflections);
		break;

	case 60:
		geometry6OCotherSymmetries(allReflections);
		break;

	case 61:
		geometry6TPRotherSymmetries(allReflections);
		break;

	case 62:
		geometry6HPotherSymmetries(allReflections);
		break;

	case 63:
		geometry6PPYotherSymmetries(allReflections);
		break;

	case 70:
		geometry7COCotherSymmetries(allReflections);
		break;

	case 71:
		geometry7PBPYotherSymmetries(allReflections);
		break;

	case 72:
		geometry7CTPRotherSymmetries(allReflections);
		break;

	case 73:
		geometry7HPYotherSymmetries(allReflections);
		break;

	case 74:
		geometry7HPotherSymmetries(allReflections);
		break;

	case 75:
		geometry7JETPYotherSymmetries(allReflections);
		break;

	case 80:
		geometry8SAPRotherSymmetries(allReflections);
		break;

	case 81:
		geometry8TDDotherSymmetries(allReflections);
		break;

	case 82:
		geometry8BTPRotherSymmetries(allReflections);
		break;

	case 83:
		geometry8HBPYotherSymmetries(allReflections);
		break;

	case 84:
		geometry8CUotherSymmetries(allReflections);
		break;

	case 85:
		geometry8ETBPYotherSymmetries(allReflections);
		break;

	case 86:
		geometry8HPYotherSymmetries(allReflections);
		break;

	case 87:
		geometry8OPotherSymmetries(allReflections);
		break;

	case 88:
		geometry8JGBFotherSymmetries(allReflections);
		break;

	case 90:
		geometry9TCTPRotherSymmetries(allReflections);
		break;

	case 91:
		geometry9CSAPRotherSymmetries(allReflections);
		break;

	case 92:
		geometry9MFFotherSymmetries(allReflections);
		break;

	case 93:
		geometry9CCUotherSymmetries(allReflections);
		break;

	case 94:
		geometry9HHotherSymmetries(allReflections);
		break;

	case 95:
		geometry9OPYotherSymmetries(allReflections);
		break;

	case 96:
		geometry9EPotherSymmetries(allReflections);
		break;

	case 97:
		geometry9HBPYotherSymmetries(allReflections);
		break;

	case 98:
		geometry9JTCotherSymmetries(allReflections);
		break;

	case 99:
		geometry9JTDICotherSymmetries(allReflections);
		break;

	case 100:
		//return geometry10PointSphere(mol0, cutAngle, reflectionOperation);
		break;

	case 101:
		//return geometry10TD(mol0, cutAngle, reflectionOperation);
		break;

	case 102:
		//return geometry10JSPC(mol0, cutAngle, reflectionOperation);
		break;

	case 103:
		//return geometry10JBCSAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 110:
		//return geometry11JCPAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 120:
		//return geometry12IC(mol0, cutAngle, reflectionOperation);
		break;

	default:
		cout << "Geometry not found" << endl;
		exit(1);
		break;
	}
}


std::string Geometries::selectGeometrySymmetriesFlag(
	int select,
	int iSymmetry,
	int symmetryType)
{
	switch (select)
	{
	case 40:
		return geometry4TetrahedronSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 41:
		return geometry4SquareSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 42:
		return geometry4SSSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 43:
		return geometry4vTBPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 50:
		return geometry5TBPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 51:
		return geometry5SPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 52:
		//return geometry5VOC(mol0, cutAngle, reflectionOperation);
		break;

	case 53:
		return geometry5PPSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 60:
		return geometry6OCSymmetryFlags(iSymmetry,symmetryType);
		break;

	case 61:
		return geometry6TPRSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 62:
		return geometry6HPSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 63:
		return geometry6PPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 70:
		return geometry7COCSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 71:
		return geometry7PBPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 72:
		return geometry7CTPRSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 73:
		return geometry7HPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 74:
		return geometry7HPSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 75:
		return geometry7JETPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 80:
		return geometry8SAPRSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 81:
		return geometry8TDDSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 82:
		return geometry8BTPRSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 83:
		return geometry8HBPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 84:
		return geometry8CUSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 85:
		return geometry8ETBPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 86:
		return geometry8HPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 87:
		return geometry8OPSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 88:
		return geometry8JGBFSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 90:
		return geometry9TCTPRSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 91:
		return geometry9CSAPRSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 92:
		return geometry9MFFSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 93:
		return geometry9CCUSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 94:
		return geometry9HHSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 95:
		return geometry9OPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 96:
		return geometry9EPSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 97:
		return geometry9HBPYSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 98:
		return geometry9JTCSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 99:
		return geometry9JTDICSymmetryFlags(iSymmetry, symmetryType);
		break;

	case 100:
		//return geometry10PointSphere(mol0, cutAngle, reflectionOperation);
		break;

	case 101:
		cout << "OTHER SYMMETRIES NOT DONE " << endl;
		exit(1);
		//return geometry10TD(mol0, cutAngle, reflectionOperation);
		break;

	case 102:
		//return geometry10JSPC(mol0, cutAngle, reflectionOperation);
		break;

	case 103:
		//return geometry10JBCSAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 110:
		//return geometry11JCPAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 120:
		//return geometry12IC(mol0, cutAngle, reflectionOperation);
		break;

	default:
		cout << "Geometry not found" << endl;
		exit(1);
		break;
	}
}



string Geometries::sizeToGeometryCode(int geoCode)
{
	switch (geoCode)
	{
	case 40:
		return "T-4";
		break;

	case 41:
		return "SP-4";
		break;

	case 42:
		return "SS-4";
		break;

	case 43:
		return "vTBPY-4";
		break;

	case 50:
		return "TBPY-5";
		break;

	case 51:
		return "SPY-5";
		break;

	case 52:
		return "vOC-5";
		break;

	case 53:
		return "PP-5";
		break;

	case 60:
		return "OC-6";
		break;

	case 61:
		return "TPR-6";
		break;

	case 62:
		return "HP-6";
		break;

	case 63:
		return "PPY-6";
		break;

	case 70:
		return "COC-7";
		break;

	case 71:
		return "PBPY-7";
		break;

	case 72:
		return "CTPR-7";
		break;

	case 73:
		return "HPY-7";
		break;

	case 74:
		return "HP-7";
		break;

	case 75:
		return "JETPY-7";
		break;

	case 80:
		return "SAPR-8";
		break;

	case 81:
		return "TDD-8";
		break;

	case 82:
		return "BTPR-8";
		break;

	case 83:
		return "HBPY-8";
		break;

	case 84:
		return "CU-8";
		break;

	case 85:
		return "ETBPY-8";
		break;

	case 86:
		return "HPY-8";
		break;

	case 87:
		return "OP-8";
		break;

	case 88:
		return "JGBF-8";
		break;

	case 90:
		return "TCTPR-9";
		break;

	case 91:
		return "CSAPR-9";
		break;

	case 92:
		return "MFF-9";
		break;

	case 93:
		return "CCU-9";
		break;

	case 94:
		return "HH-9";
		break;

	case 95:
		return "OPY-9";
		break;

	case 96:
		return "EP-9";
		break;

	case 97:
		return "HBPY-9";
		break;

	case 98:
		return "JTC-9";
		break;

	case 99:
		return "JTDIC-9";
		break;

	case 100:
		return "PointSphere-10";
		break;

	case 101:
		return "TD-10";
		break;

	case 102:
		return "JSPC-10";
		break;

	case 103:
		return "JBCSAPR-10";
		break;

	case 110:
		return "JCPAPR-11";
		break;

	case 120:
		return "IC-12";
		break;

	default:
		cout << "Geometry not found" << endl;
		exit(1);
		break;
	}


}


std::vector<double> Geometries::geometry4Tetrahedron(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 4;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;
	vector<double> vectorRotations(44);

	mol0[0].x = 0.000000;
	mol0[0].y = 0.816300;
	mol0[0].z = -0.577200;
	mol0[1].x = 0.000000;
	mol0[1].y = -0.816300;
	mol0[1].z = -0.577200;
	mol0[2].x = 0.816300;
	mol0[2].y = 0.000000;
	mol0[2].z = 0.577200;
	mol0[3].x = -0.816300;
	mol0[3].y = 0.000000;
	mol0[3].z = 0.577200;

	vector<double> auxReferenceAxis(3);
	// ROTATIONS
	//C3-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = -2.0e0 * auxMath_._pi / 3.0e0;
	//C3-2-p
	vectorRotations[8] = mol0[1].x;
	vectorRotations[9] = mol0[1].y;
	vectorRotations[10] = mol0[1].z;
	vectorRotations[11] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-2-m
	vectorRotations[12] = mol0[1].x;
	vectorRotations[13] = mol0[1].y;
	vectorRotations[14] = mol0[1].z;
	vectorRotations[15] = -2.0e0 * auxMath_._pi / 3.0e0;
	//C3-3-p
	vectorRotations[16] = mol0[2].x;
	vectorRotations[17] = mol0[2].y;
	vectorRotations[18] = mol0[2].z;
	vectorRotations[19] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-3-m
	vectorRotations[20] = mol0[2].x;
	vectorRotations[21] = mol0[2].y;
	vectorRotations[22] = mol0[2].z;
	vectorRotations[23] = -2.0e0 * auxMath_._pi / 3.0e0;
	//C3-4-p
	vectorRotations[24] = mol0[3].x;
	vectorRotations[25] = mol0[3].y;
	vectorRotations[26] = mol0[3].z;
	vectorRotations[27] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-4-m
	vectorRotations[28] = mol0[3].x;
	vectorRotations[29] = mol0[3].y;
	vectorRotations[30] = mol0[3].z;
	vectorRotations[31] = -2.0e0 * auxMath_._pi / 3.0e0;

	//C2-1
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = auxMath_._pi;
	//C2-2
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[36] = auxReferenceAxis[0];
	vectorRotations[37] = auxReferenceAxis[1];
	vectorRotations[38] = auxReferenceAxis[2];
	vectorRotations[39] = auxMath_._pi;
	//C3-3
	auxReferenceAxis[0] = 0.5e0 * (mol0[3].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[3].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[3].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[40] = auxReferenceAxis[0];
	vectorRotations[41] = auxReferenceAxis[1];
	vectorRotations[42] = auxReferenceAxis[2];
	vectorRotations[43] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[2] = 3;
	reflectionOperation[3] = 2;

	return vectorRotations;
}


std::vector<double> Geometries::geometry4Square(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 4;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;
	vector<double> vectorRotations(28);

	mol0[0].x = 1.000000000000000;
	mol0[0].y = 0.000000000000000;
	mol0[0].z = 0.000000000000000;
	mol0[1].x = 0.000000000000000;
	mol0[1].y = 1.000000000000000;
	mol0[1].z = 0.000000000000000;
	mol0[2].x = -1.000000000000000;
	mol0[2].y = 0.000000000000000;
	mol0[2].z = 0.000000000000000;
	mol0[3].x = 0.000000000000000;
	mol0[3].y = -1.000000000000000;
	mol0[3].z = 0.000000000000000;

	vector<double> auxReferenceAxis(3);
	// ROTATIONS
	//C2-2 ( C4-1-p C2-1 C4-1-m C2-4 )
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi;
	//C2-3 ( C4-1-p C2-1 C4-1-m C2-5 )
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	//C3-4 ( C4-1-p C2-1 C4-1-m C2-2 )
	vectorRotations[8] = mol0[1].x;
	vectorRotations[9] = mol0[1].y;
	vectorRotations[10] = mol0[1].z;
	vectorRotations[11] = auxMath_._pi;
	//C2-5 ( C4-1-p C2-1 C4-1-m C2-3 )
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;
	// C4-1-p
	vectorRotations[16] = 0.00000000000000000;
	vectorRotations[17] = 0.00000000000000000;
	vectorRotations[18] = 1.00000000000000000;
	vectorRotations[19] = auxMath_._pi / 2.0e0;
	// C2-1 ( C2-2 C2-3 C2-4 C2-5 )
	vectorRotations[20] = 0.00000000000000000;
	vectorRotations[21] = 0.00000000000000000;
	vectorRotations[22] = 1.00000000000000000;
	vectorRotations[23] = auxMath_._pi;
	// C4-1-m
	vectorRotations[24] = 0.00000000000000000;
	vectorRotations[25] = 0.00000000000000000;
	vectorRotations[26] = 1.00000000000000000;
	vectorRotations[27] = 3.0e0 * auxMath_._pi / 2.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 2;
	reflectionOperation[2] = 0;

	return vectorRotations;
}

std::vector<double> Geometries::geometry4SS(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 4;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(4);

	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.00000000;
	mol0[2].y = 1.00000000;
	mol0[2].z = 0.00000000;
	mol0[3].x = 0.00000000;
	mol0[3].y = 0.00000000;
	mol0[3].z = 1.00000000;

	//C2-1 
	vector<double> auxReferenceAxis(3);
	auxReferenceAxis[0] = mol0[1].x + mol0[2].x;
	auxReferenceAxis[1] = mol0[1].y + mol0[2].y;
	auxReferenceAxis[2] = mol0[1].z + mol0[2].z;
	auxMath_.normalize(auxReferenceAxis);

	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 2;
	reflectionOperation[2] = 1;

	return vectorRotations;
}

std::vector<double> Geometries::geometry4vTBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 4;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(8);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = -0.50000000;
	mol0[2].y = 0.86602540;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.50000000;
	mol0[3].y = -0.86602540;
	mol0[3].z = 0.00000000;

	//C3-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 3.0e0);
	//C3-1-m
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 3.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 3;
	reflectionOperation[3] = 1;

	return vectorRotations;
}


std::vector<double> Geometries::geometry5TBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 5;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 1.9e0;
	vector<double> vectorRotations(20);

	mol0[0].x = 0.0000000000;
	mol0[0].y = 0.0000000000;
	mol0[0].z = 1.0000000000;
	mol0[1].x = 1.0000000000;
	mol0[1].y = 0.0000000000;
	mol0[1].z = 0.0000000000;
	mol0[2].x = -1.000000000;
	mol0[2].y = 0.0000000000;
	mol0[2].z = 0.0000000000;
	mol0[3].x = 0.0000000000;
	mol0[3].y = 0.8660254000;
	mol0[3].z = -0.500000000;
	mol0[4].x = 0.0000000000;
	mol0[4].y = -0.866025400;
	mol0[4].z = -0.500000000;

	//C3-1-p
	vectorRotations[0] = 1.0e0;
	vectorRotations[1] = 0.0e0;
	vectorRotations[2] = 0.0e0;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m
	vectorRotations[4] = 1.0e0;
	vectorRotations[5] = 0.0e0;
	vectorRotations[6] = 0.0e0;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;
	//C2-1 ( C3-1-p C3-1-m )
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = auxMath_._pi;
	//C2-2 ( C3-1-p C3-1-m )
	vectorRotations[12] = mol0[3].x;
	vectorRotations[13] = mol0[3].y;
	vectorRotations[14] = mol0[3].z;
	vectorRotations[15] = auxMath_._pi;
	//C2-3 ( C3-1-p C3-1-m )
	vectorRotations[16] = mol0[4].x;
	vectorRotations[17] = mol0[4].y;
	vectorRotations[18] = mol0[4].z;
	vectorRotations[19] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 2;
	reflectionOperation[2] = 1;

	return vectorRotations;
}


std::vector<double> Geometries::geometry5SPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 5;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;
	vector<double> vectorRotations(12);


	mol0[0].x = 0.0000000000;
	mol0[0].y = 0.0000000000;
	mol0[0].z = 1.0000000000;
	mol0[1].x = 0.968245837;
	mol0[1].y = 0.0000000000;
	mol0[1].z = -0.250000000;
	mol0[2].x = -0.968245837;
	mol0[2].y = 0.0000000000;
	mol0[2].z = -0.250000000;
	mol0[3].x = 0.0000000000;
	mol0[3].y = 0.9682458370;
	mol0[3].z = -0.250000000;
	mol0[4].x = 0.0000000000;
	mol0[4].y = -0.968245837;
	mol0[4].z = -0.250000000;

	//C4-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C2-1 
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = auxMath_._pi;
	//C4-1-m
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = -auxMath_._pi / 2.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 2;
	reflectionOperation[2] = 1;

	return vectorRotations;
}


std::vector<double> Geometries::geometry5VOC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 5;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 0.9 * auxMath_._pi;
	vector<double> vectorRotations(12);


	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.00000000;
	mol0[2].y = 1.00000000;
	mol0[2].z = 0.00000000;
	mol0[3].x = -1.00000000;
	mol0[3].y = 0.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = 0.00000000;
	mol0[4].y = -1.00000000;
	mol0[4].z = 0.00000000;

	//c4 - 1
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//c4 - 1
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = auxMath_._pi;
	//c4 - 1
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = -auxMath_._pi / 2.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 3;
	reflectionOperation[3] = 1;

	return vectorRotations;
}

std::vector<double> Geometries::geometry5PP(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 5;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(36);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.30901699;	
	mol0[1].y = 0.95105652;
	mol0[1].z = 0.00000000;
	mol0[2].x = -0.80901699;
	mol0[2].y = 0.58778525;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.80901699;
	mol0[3].y = -0.58778525;
	mol0[3].z = 0.00000000;
	mol0[4].x = 0.30901699;
	mol0[4].y = -0.95105652;
	mol0[4].z = 0.00000000;

	auxReferenceAxis[0] = 0.0e0;
	auxReferenceAxis[1] = 0.0e0;
	auxReferenceAxis[2] = 1.0e0;

	//C5-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 5.0e0);
	//C5-1-pp
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
	//C5-1-mm
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
	//C5-1-m
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);

	//C2-1 ( C5-1-p C5-1-pp C5-1-mm C5-1-m )
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = auxMath_._pi;
	//C2-2 ( C5-1-p C5-1-pp C5-1-mm C5-1-m )
	vectorRotations[20] = mol0[1].x;
	vectorRotations[21] = mol0[1].y;
	vectorRotations[22] = mol0[1].z;
	vectorRotations[23] = auxMath_._pi;
	//C2-3 ( C5-1-p C5-1-pp C5-1-mm C5-1-m )
	vectorRotations[24] = mol0[2].x;
	vectorRotations[25] = mol0[2].y;
	vectorRotations[26] = mol0[2].z;
	vectorRotations[27] = auxMath_._pi;
	//C2-4 ( C5-1-p C5-1-pp C5-1-mm C5-1-m )
	vectorRotations[28] = mol0[3].x;
	vectorRotations[29] = mol0[3].y;
	vectorRotations[30] = mol0[3].z;
	vectorRotations[31] = auxMath_._pi;
	//C2-5 ( C5-1-p C5-1-pp C5-1-mm C5-1-m )
	vectorRotations[32] = mol0[4].x;
	vectorRotations[33] = mol0[4].y;
	vectorRotations[34] = mol0[4].z;
	vectorRotations[35] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;

	return vectorRotations;
}


std::vector<double> Geometries::geometry6OC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 6;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;
	vector<double> vectorRotations(92);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.00000000;
	mol0[2].y = 1.00000000;
	mol0[2].z = 0.00000000;
	mol0[3].x = -1.00000000;
	mol0[3].y = 0.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = 0.00000000;
	mol0[4].y = -1.00000000;
	mol0[4].z = 0.00000000;
	mol0[5].x = 0.00000000;
	mol0[5].y = 0.00000000;
	mol0[5].z = -1.00000000;

	//C4-1-p [ref 0]
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C2-1 ( C4-2-p C4-2-m C2-2 C4-3-p C4-3-m C2-3 C2-8 C2-9 ) [ref 0]
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = auxMath_._pi;
	//C4-1-m [ref 0]
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * auxMath_._pi / 2.0e0;

	//C4-2-p  [ref 1]
	vectorRotations[12] = mol0[1].x;
	vectorRotations[13] = mol0[1].y;
	vectorRotations[14] = mol0[1].z;
	vectorRotations[15] = auxMath_._pi / 2.0e0;
	//C2-2 ( C4-1-p C4-1-m C2-1 C4-3-p C4-3-m C2-3 C2-5 C2-7 )  [ref 1]
	vectorRotations[16] = mol0[1].x;
	vectorRotations[17] = mol0[1].y;
	vectorRotations[18] = mol0[1].z;
	vectorRotations[19] = auxMath_._pi;
	//C4-2-m  [ref 1]
	vectorRotations[20] = mol0[1].x;
	vectorRotations[21] = mol0[1].y;
	vectorRotations[22] = mol0[1].z;
	vectorRotations[23] = 3.0e0 * auxMath_._pi / 2.0e0;

	//C4-3-p  [ref 2]
	vectorRotations[24] = mol0[2].x;
	vectorRotations[25] = mol0[2].y;
	vectorRotations[26] = mol0[2].z;
	vectorRotations[27] = auxMath_._pi / 2.0e0;
	//C2-3 ( C4-1-p C4-1-m C2-1 C4-2-p C4-2-m C2-2 C2-4 C2-6 )  [ref 2]
	vectorRotations[28] = mol0[2].x;
	vectorRotations[29] = mol0[2].y;
	vectorRotations[30] = mol0[2].z;
	vectorRotations[31] = auxMath_._pi;
	//C4-3-m  [ref 2]
	vectorRotations[32] = mol0[2].x;
	vectorRotations[33] = mol0[2].y;
	vectorRotations[34] = mol0[2].z;
	vectorRotations[35] = 3.0e0 * auxMath_._pi / 2.0e0;

	//C2-4 ( C4-3-p C4-3-m C2-3 C3-2-p C3-2-m C3-3-p C3-3-m C2-6 )  [0-1]
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[36] = auxReferenceAxis[0];
	vectorRotations[37] = auxReferenceAxis[1];
	vectorRotations[38] = auxReferenceAxis[2];
	vectorRotations[39] = auxMath_._pi;

	//C2-5 ( C4-2-p C4-2-m C2-2 C3-3-p C3-3-m C3-4-p C3-4-m C2-7 )  [0-2]
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[40] = auxReferenceAxis[0];
	vectorRotations[41] = auxReferenceAxis[1];
	vectorRotations[42] = auxReferenceAxis[2];
	vectorRotations[43] = auxMath_._pi;

	//C2-6 ( C4-3-p C4-3-m C2-3 C3-1-p C3-1-m C3-4-p C3-4-m C2-4 )[0-3]
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[44] = auxReferenceAxis[0];
	vectorRotations[45] = auxReferenceAxis[1];
	vectorRotations[46] = auxReferenceAxis[2];
	vectorRotations[47] = auxMath_._pi;

	//C2-7  ( C4-2-p C4-2-m C2-2 C3-1-p C3-1-m C3-2-p C3-2-m C2-5 )  [0-4]
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[48] = auxReferenceAxis[0];
	vectorRotations[49] = auxReferenceAxis[1];
	vectorRotations[50] = auxReferenceAxis[2];
	vectorRotations[51] = auxMath_._pi;

	//C2-8  ( C4-1-p C4-1-m C2-1 C3-2-p C3-2-m C3-4-p C3-4-m C2-9 ) [1-2]
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[52] = auxReferenceAxis[0];
	vectorRotations[53] = auxReferenceAxis[1];
	vectorRotations[54] = auxReferenceAxis[2];
	vectorRotations[55] = auxMath_._pi;

	//C2-9  ( C4-1-p C4-1-m C2-1 C3-1-p C3-1-m C3-3-p C3-3-m C2-8 )  [1-4]
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[56] = auxReferenceAxis[0];
	vectorRotations[57] = auxReferenceAxis[1];
	vectorRotations[58] = auxReferenceAxis[2];
	vectorRotations[59] = auxMath_._pi;

	//C3-1-p [0-1-2]
	auxReferenceAxis[0] = (mol0[0].x + mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[60] = auxReferenceAxis[0];
	vectorRotations[61] = auxReferenceAxis[1];
	vectorRotations[62] = auxReferenceAxis[2];
	vectorRotations[63] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m [0-1-2]
	vectorRotations[64] = auxReferenceAxis[0];
	vectorRotations[65] = auxReferenceAxis[1];
	vectorRotations[66] = auxReferenceAxis[2];
	vectorRotations[67] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C3-2-p  [0-2-3]
	auxReferenceAxis[0] = (mol0[0].x + mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[68] = auxReferenceAxis[0];
	vectorRotations[69] = auxReferenceAxis[1];
	vectorRotations[70] = auxReferenceAxis[2];
	vectorRotations[71] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-2-m  [0-2-3]
	vectorRotations[72] = auxReferenceAxis[0];
	vectorRotations[73] = auxReferenceAxis[1];
	vectorRotations[74] = auxReferenceAxis[2];
	vectorRotations[75] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C3-3-p  [0-3-4]
	auxReferenceAxis[0] = (mol0[0].x + mol0[3].x + mol0[4].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[3].y + mol0[4].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[3].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[76] = auxReferenceAxis[0];
	vectorRotations[77] = auxReferenceAxis[1];
	vectorRotations[78] = auxReferenceAxis[2];
	vectorRotations[79] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-3-m  [0-3-4]
	vectorRotations[80] = auxReferenceAxis[0];
	vectorRotations[81] = auxReferenceAxis[1];
	vectorRotations[82] = auxReferenceAxis[2];
	vectorRotations[83] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C3-4-p  [0-1-4]
	auxReferenceAxis[0] = (mol0[0].x + mol0[4].x + mol0[1].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[4].y + mol0[1].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[4].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[84] = auxReferenceAxis[0];
	vectorRotations[85] = auxReferenceAxis[1];
	vectorRotations[86] = auxReferenceAxis[2];
	vectorRotations[87] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-4-m  [0-1-4]
	vectorRotations[88] = auxReferenceAxis[0];
	vectorRotations[89] = auxReferenceAxis[1];
	vectorRotations[90] = auxReferenceAxis[2];
	vectorRotations[91] = 4.0e0 * auxMath_._pi / 3.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 5;
	reflectionOperation[5] = 0;

	return vectorRotations;
}


std::vector<double> Geometries::geometry6TPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 6;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;
	vector<double> vectorRotations(20);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = -0.65465367070798;
	mol0[0].y = -0.37796447300923;
	mol0[0].z = 0.65465367070798;
	mol0[1].x = -0.65465367070798;
	mol0[1].y = -0.37796447300923;
	mol0[1].z = -0.65465367070798;
	mol0[2].x = 0.65465367070798;
	mol0[2].y = -0.37796447300923;
	mol0[2].z = 0.65465367070798;
	mol0[3].x = 0.65465367070798;
	mol0[3].y = -0.37796447300923;
	mol0[3].z = -0.65465367070798;
	mol0[4].x = 0.00000000000000;
	mol0[4].y = 0.75592894601846;
	mol0[4].z = 0.65465367070798;
	mol0[5].x = 0.00000000000000;
	mol0[5].y = 0.75592894601846;
	mol0[5].z = -0.65465367070798;

	//C3-1-p [ref 0]
	vectorRotations[0] = 0.000000000;
	vectorRotations[1] = 0.000000000;
	vectorRotations[2] = 1.000000000;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m [ref 0]
	vectorRotations[4] = 0.000000000;
	vectorRotations[5] = 0.000000000;
	vectorRotations[6] = 1.000000000;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C2-1
	auxReferenceAxis[0] = 0.5e0 * (mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = auxMath_._pi;

	//C2-2
	auxReferenceAxis[0] = 0.5e0 * (mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;

	//C2-3
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 4;
	reflectionOperation[4] = 0;
	reflectionOperation[1] = 5;
	reflectionOperation[5] = 1;

	return vectorRotations;
}

std::vector<double> Geometries::geometry6HP(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 6;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(44);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 1.00000000;	
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.50000000;
	mol0[1].y = 0.86602540;
	mol0[1].z = 0.00000000;
	mol0[2].x = -0.50000000;
	mol0[2].y = 0.86602540;
	mol0[2].z = 0.00000000;
	mol0[3].x = -1.00000000;
	mol0[3].y = 0.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.50000000;
	mol0[4].y = -0.86602540;
	mol0[4].z = 0.00000000;
	mol0[5].x = 0.50000000;
	mol0[5].y = -0.86602540;
	mol0[5].z = 0.00000000;

	auxReferenceAxis[0] = 0.0e0;
	auxReferenceAxis[1] = 0.0e0;
	auxReferenceAxis[2] = 1.0e0;

	//C6-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 6.0e0);
	//C3-1-p
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);
	//C2-1 ( C2-2 C2-3 C2-4 C2-5 C2-6 C2-7 ) 
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);
	//C3-1-m
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);
	//C6-1-m
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = 5.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);



	//C2-2 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-5 )
	vectorRotations[20] = mol0[0].x;
	vectorRotations[21] = mol0[0].y;
	vectorRotations[22] = mol0[0].z;
	vectorRotations[23] = auxMath_._pi;

	//C2-3 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-6 )
	auxReferenceAxis[0] = mol0[0].x + mol0[1].x;
	auxReferenceAxis[1] = mol0[0].y + mol0[1].y;
	auxReferenceAxis[2] = mol0[0].z + mol0[1].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = auxMath_._pi;

	//C2-4 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-7 )
	vectorRotations[28] = mol0[1].x;
	vectorRotations[29] = mol0[1].y;
	vectorRotations[30] = mol0[1].z;
	vectorRotations[31] = auxMath_._pi;

	//C2-5 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-2 )
	auxReferenceAxis[0] = mol0[1].x + mol0[2].x;
	auxReferenceAxis[1] = mol0[1].y + mol0[2].y;
	auxReferenceAxis[2] = mol0[1].z + mol0[2].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = auxMath_._pi;

	//C2-6 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-3 )
	vectorRotations[36] = mol0[2].x;
	vectorRotations[37] = mol0[2].y;
	vectorRotations[38] = mol0[2].z;
	vectorRotations[39] = auxMath_._pi;

	//C2-7 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-4 )
	auxReferenceAxis[0] = mol0[2].x + mol0[3].x;
	auxReferenceAxis[1] = mol0[2].y + mol0[3].y;
	auxReferenceAxis[2] = mol0[2].z + mol0[3].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[40] = auxReferenceAxis[0];
	vectorRotations[41] = auxReferenceAxis[1];
	vectorRotations[42] = auxReferenceAxis[2];
	vectorRotations[43] = auxMath_._pi;


	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 4;
	reflectionOperation[4] = 0;
	reflectionOperation[1] = 5;
	reflectionOperation[5] = 1;

	return vectorRotations;
}


std::vector<double> Geometries::geometry6PPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 6;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(16);
	vector<double> auxReferenceAxis(3);
	
	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.30901699;
	mol0[2].y = 0.95105652;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.80901699;
	mol0[3].y = 0.58778525;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.80901699;
	mol0[4].y = -0.58778525;
	mol0[4].z = 0.00000000;
	mol0[5].x = 0.30901699;
	mol0[5].y = -0.95105652;
	mol0[5].z = 0.00000000;

	//C5-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 5.0e0);
	//C5-1-pp
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
	//C5-1-mm
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
	//C5-1-m
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[3] = 4;
	reflectionOperation[4] = 3;
	reflectionOperation[2] = 5;
	reflectionOperation[5] = 2;

	return vectorRotations;
}


std::vector<double> Geometries::geometry7COC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;
	vector<double> vectorRotations(8);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 1.00000000;
	mol0[1].x = 0.97767167;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.21013831;
	mol0[2].x = 0.16977090;
	mol0[2].y = 0.96281864;
	mol0[2].z = 0.21013831;
	mol0[3].x = -0.91871085;
	mol0[3].y = 0.33438340;
	mol0[3].z = 0.21013831;
	mol0[4].x = -0.48883583;
	mol0[4].y = -0.84668850;
	mol0[4].z = 0.21013831;
	mol0[5].x = 0.36282725;
	mol0[5].y = -0.62843523;
	mol0[5].z = -0.68805926;
	mol0[6].x = -0.26010411;
	mol0[6].y = 0.45051354;
	mol0[6].z = -0.85403946;

	//C3-1-p
	vectorRotations[0] = mol0[5].x;
	vectorRotations[1] = mol0[5].y;
	vectorRotations[2] = mol0[5].z;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m
	vectorRotations[4] = mol0[5].x;
	vectorRotations[5] = mol0[5].y;
	vectorRotations[6] = mol0[5].z;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 2;
	reflectionOperation[2] = 0;
	reflectionOperation[4] = 6;
	reflectionOperation[6] = 4;


	return vectorRotations;
}

std::vector<double> Geometries::geometry7PBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;;
	vector<double> vectorRotations(36);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.000000000000;
	mol0[0].y = 0.000000000000;
	mol0[0].z = 1.000000000000;
	mol0[1].x = 0.000000000000;
	mol0[1].y = 0.000000000000;
	mol0[1].z = -1.000000000000;
	mol0[2].x = 1.000000000000;
	mol0[2].y = 0.000000000000;
	mol0[2].z = 0.000000000000;
	mol0[3].x = 0.309016994375;
	mol0[3].y = 0.951056516295;
	mol0[3].z = 0.000000000000;
	mol0[4].x = -0.809016994375;
	mol0[4].y = 0.587785252292;
	mol0[4].z = 0.000000000000;
	mol0[5].x = -0.809016994375;
	mol0[5].y = -0.587785252292;
	mol0[5].z = 0.000000000000;
	mol0[6].x = 0.309016994375;
	mol0[6].y = -0.951056516295;
	mol0[6].z = 0.000000000000;

	//C5-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 5.0e0;
	//C5-1-pp
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 5.0e0;
	//C5-1-m
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 6.0e0 * auxMath_._pi / 5.0e0;
	//C5-1-mm
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 8.0e0 * auxMath_._pi / 5.0e0;

	//C2-1 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )
	vectorRotations[16] = mol0[2].x;
	vectorRotations[17] = mol0[2].y;
	vectorRotations[18] = mol0[2].z;
	vectorRotations[19] = auxMath_._pi;
	//C2-2 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )
	vectorRotations[20] = mol0[3].x;
	vectorRotations[21] = mol0[3].y;
	vectorRotations[22] = mol0[3].z;
	vectorRotations[23] = auxMath_._pi;
	//C2-3 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )
	vectorRotations[24] = mol0[4].x;
	vectorRotations[25] = mol0[4].y;
	vectorRotations[26] = mol0[4].z;
	vectorRotations[27] = auxMath_._pi;
	//C2-4 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )
	vectorRotations[28] = mol0[5].x;
	vectorRotations[29] = mol0[5].y;
	vectorRotations[30] = mol0[5].z;
	vectorRotations[31] = auxMath_._pi;
	//C2-5 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )
	vectorRotations[32] = mol0[6].x;
	vectorRotations[33] = mol0[6].y;
	vectorRotations[34] = mol0[6].z;
	vectorRotations[35] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 1;
	reflectionOperation[1] = 0;

	return vectorRotations;
}

std::vector<double> Geometries::geometry7CTPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;;
	vector<double> vectorRotations(4);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.0000;
	mol0[0].y = 0.0000;
	mol0[0].z = 1.0000;
	mol0[1].x = 0.6869;
	mol0[1].y = 0.6869;
	mol0[1].z = 0.2374;
	mol0[2].x = -0.6869;
	mol0[2].y = 0.6869;
	mol0[2].z = 0.2374;
	mol0[3].x = 0.6869;
	mol0[3].y = -0.6869;
	mol0[3].z = 0.2374;
	mol0[4].x = -0.6869;
	mol0[4].y = -0.6869;
	mol0[4].z = 0.2374;
	mol0[5].x = 0.6175;
	mol0[5].y = 0.0000;
	mol0[5].z = -0.7866;
	mol0[6].x = -0.6175;
	mol0[6].y = 0.0000;
	mol0[6].z = -0.7866;

	//C2-1 
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 3;
	reflectionOperation[3] = 1;
	reflectionOperation[2] = 4;
	reflectionOperation[4] = 2;

	return vectorRotations;
}

std::vector<double> Geometries::geometry7HPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(20);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;	
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.50000000;
	mol0[2].y = 0.86602540;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.50000000;
	mol0[3].y = 0.86602540;
	mol0[3].z = 0.00000000;
	mol0[4].x = -1.00000000;
	mol0[4].y = 0.00000000;
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.50000000;
	mol0[5].y = -0.86602540;
	mol0[5].z = 0.00000000;
	mol0[6].x = 0.50000000;
	mol0[6].y = -0.86602540;
	mol0[6].z = 0.00000000;

	//C6-1-p 
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 6.0e0);
	//C3-1-p 
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);
	//C2-1 
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);
	//C3-1-m 
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);
	//C6-1-m 
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = 5.0e0 * (2.0e0 * auxMath_._pi / 6.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[2] = 6;
	reflectionOperation[6] = 2;
	reflectionOperation[3] = 5;
	reflectionOperation[5] = 3;

	return vectorRotations;
}


std::vector<double> Geometries::geometry7HP(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(52);
	vector<double> auxReferenceAxis(3);	

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;	
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.62348980;
	mol0[1].y = 0.78183148;
	mol0[1].z = 0.00000000;
	mol0[2].x = -0.22252093;
	mol0[2].y = 0.97492791;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.90096887;
	mol0[3].y = 0.43388374;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.90096887;
	mol0[4].y = -0.43388374;
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.22252093;
	mol0[5].y = -0.97492791;
	mol0[5].z = 0.00000000;
	mol0[6].x = 0.62348980;
	mol0[6].y = -0.78183148;
	mol0[6].z = 0.00000000;

	auxReferenceAxis[0] = 0.0e0;
	auxReferenceAxis[1] = 0.0e0;
	auxReferenceAxis[2] = -2.54417e0;
	auxMath_.normalize(auxReferenceAxis);

	//C7-1-p 
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-pp 
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-ppp 
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-mmm 
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-mm 
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = 5.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-m 
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = 6.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);

	//C2-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[24] = mol0[0].x;
	vectorRotations[25] = mol0[0].y;
	vectorRotations[26] = mol0[0].z;
	vectorRotations[27] = auxMath_._pi;
	//C2-2 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[28] = mol0[1].x;
	vectorRotations[29] = mol0[1].y;
	vectorRotations[30] = mol0[1].z;
	vectorRotations[31] = auxMath_._pi;
	//C2-3 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[32] = mol0[2].x;
	vectorRotations[33] = mol0[2].y;
	vectorRotations[34] = mol0[2].z;
	vectorRotations[35] = auxMath_._pi;
	//C2-4 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[36] = mol0[3].x;
	vectorRotations[37] = mol0[3].y;
	vectorRotations[38] = mol0[3].z;
	vectorRotations[39] = auxMath_._pi;
	//C2-5 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[40] = mol0[4].x;
	vectorRotations[41] = mol0[4].y;
	vectorRotations[42] = mol0[4].z;
	vectorRotations[43] = auxMath_._pi;
	//C2-6 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[44] = mol0[5].x;
	vectorRotations[45] = mol0[5].y;
	vectorRotations[46] = mol0[5].z;
	vectorRotations[47] = auxMath_._pi;
	//C2-7 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m )
	vectorRotations[48] = mol0[6].x;
	vectorRotations[49] = mol0[6].y;
	vectorRotations[50] = mol0[6].z;
	vectorRotations[51] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;

	return vectorRotations;
}


std::vector<double> Geometries::geometry7JETPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(8);
	vector<double> auxReferenceAxis(3);
	
	mol0[0].x = 0.75592895;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.65465367;
	mol0[1].x = 0.75592895;
	mol0[1].y = 0.00000000;
	mol0[1].z = -0.65465367;
	mol0[2].x = -0.37796447;
	mol0[2].y = 0.65465367;
	mol0[2].z = 0.65465367;
	mol0[3].x = -0.37796447;
	mol0[3].y = 0.65465367;
	mol0[3].z = -0.65465367;
	mol0[4].x = -0.37796447;
	mol0[4].y = -0.65465367;
	mol0[4].z = 0.65465367;
	mol0[5].x = -0.37796447;
	mol0[5].y = -0.65465367;
	mol0[5].z = -0.65465367;
	mol0[6].x = 0.00000000;
	mol0[6].y = 0.00000000;
	mol0[6].z = 1.00000000;

	//C3-1-p 
	vectorRotations[0] = mol0[6].x;
	vectorRotations[1] = mol0[6].y;
	vectorRotations[2] = mol0[6].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 3.0e0);
	//C3-1-m
	vectorRotations[4] = mol0[6].x;
	vectorRotations[5] = mol0[6].y;
	vectorRotations[6] = mol0[6].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 3.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 4;
	reflectionOperation[4] = 0;
	reflectionOperation[1] = 5;
	reflectionOperation[5] = 1;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8SAPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(28);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.000000;
	mol0[0].y = 0.000000;
	mol0[0].z = 1.00000000;
	mol0[1].x = 0.96528366;
	mol0[1].y = 0.000000;
	mol0[1].z = 0.26120388;
	mol0[2].x = -0.56545007;
	mol0[2].y = 0.78232905;
	mol0[2].z = 0.26120388;
	mol0[3].x = -0.88247541;
	mol0[3].y = -0.39116453;
	mol0[3].z = 0.26120387;
	mol0[4].x = 0.19991679;
	mol0[4].y = -0.94435471;
	mol0[4].z = 0.26120387;
	mol0[5].x = 0.39983358;
	mol0[5].y = 0.78232905;
	mol0[5].z = -0.47759225;
	mol0[6].x = -0.59975037;
	mol0[6].y = 0.16202565;
	mol0[6].z = -0.78361162;
	mol0[7].x = 0.48264183;
	mol0[7].y = -0.39116453;
	mol0[7].z = -0.78361162;

	// 3 4 6 7
	auxReferenceAxis[0] = 0.25e0 *
		(mol0[3].x + mol0[4].x + mol0[6].x + mol0[7].x);
	auxReferenceAxis[1] = 0.25e0 *
		(mol0[3].y + mol0[4].y + mol0[6].y + mol0[7].y);
	auxReferenceAxis[2] = 0.25e0 *
		(mol0[3].z + mol0[4].z + mol0[6].z + mol0[7].z);
	auxMath_.normalize(auxReferenceAxis);
	//C4-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C2-1 
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	//C4-1-m
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * auxMath_._pi / 2.0e0;
	//C2-2 ( C4-1-p C2-1 C4-1-m )
	auxReferenceAxis[0] = 0.5e0 *
		(mol0[1].x + mol0[7].x);
	auxReferenceAxis[1] = 0.5e0 *
		(mol0[1].y + mol0[7].y);
	auxReferenceAxis[2] = 0.5e0 *
		(mol0[1].z + mol0[7].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;
	//C2-3 ( C4-1-p C2-1 C4-1-m )
	auxReferenceAxis[0] = 0.5e0 *
		(mol0[7].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *
		(mol0[7].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *
		(mol0[7].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = auxMath_._pi;
	// C2-4 ( C4-1-p C2-1 C4-1-m )
	auxReferenceAxis[0] = 0.5e0 *
		(mol0[5].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *
		(mol0[5].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *
		(mol0[5].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = auxMath_._pi;
	//C2-5 ( C4-1-p C2-1 C4-1-m )
	auxReferenceAxis[0] = 0.5e0 *
		(mol0[6].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 *
		(mol0[6].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 *
		(mol0[6].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[2] = 5;
	reflectionOperation[5] = 2;
	reflectionOperation[1] = 0;
	reflectionOperation[0] = 1;
	reflectionOperation[7] = 3;
	reflectionOperation[3] = 7;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8TDD(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;
	vector<double> vectorRotations(12);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = -0.599725;
	mol0[0].y = 0.000000;
	mol0[0].z = 0.800200;
	mol0[1].x = -0.000025;
	mol0[1].y = -0.936400;
	mol0[1].z = 0.350900;
	mol0[2].x = 0.599775;
	mol0[2].y = 0.000000;
	mol0[2].z = 0.800200;
	mol0[3].x = -0.000025;
	mol0[3].y = 0.936400;
	mol0[3].z = 0.350900;
	mol0[4].x = -0.936425;
	mol0[4].y = 0.000000;
	mol0[4].z = -0.350900;
	mol0[5].x = -0.000025;
	mol0[5].y = -0.599700;
	mol0[5].z = -0.800200;
	mol0[6].x = 0.936475;
	mol0[6].y = 0.000000;
	mol0[6].z = -0.350900;
	mol0[7].x = -0.000025;
	mol0[7].y = 0.599700;
	mol0[7].z = -0.800200;

	//C2-1 ( C2-2 C2-3 )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;
	//C2-2 ( C2-1 C2-3 ) 
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	//C2-3 ( C2-1 C2-2 )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[7].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[7].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[7].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 2;
	reflectionOperation[2] = 0;
	reflectionOperation[4] = 6;
	reflectionOperation[6] = 4;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8BTPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(4);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = -0.654653670708;
	mol0[0].y = -0.377964473009;
	mol0[0].z = 0.654653670708;
	mol0[1].x = -0.654653670708;
	mol0[1].y = -0.377964473009;
	mol0[1].z = -0.654653670708;
	mol0[2].x = 0.654653670708;
	mol0[2].y = -0.377964473009;	
	mol0[2].z = 0.654653670708;
	mol0[3].x = 0.654653670708;
	mol0[3].y = -0.377964473009;
	mol0[3].z = -0.654653670708;
	mol0[4].x = 0.000000000000;	
	mol0[4].y = 0.755928946018;	
	mol0[4].z = 0.654653670708;
	mol0[5].x = 0.000000000000;	
	mol0[5].y = 0.755928946018;
	mol0[5].z = -0.654653670708;
	mol0[6].x = 0.000000000000;
	mol0[6].y = -1.000000000000;
	mol0[6].z = 0.000000000000;
	mol0[7].x = -0.866025403784;
	mol0[7].y = 0.500000000000;
	mol0[7].z = 0.000000000000;

	// C2-1 
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 1;
	reflectionOperation[1] = 0;
	reflectionOperation[2] = 3;
	reflectionOperation[3] = 2;
	reflectionOperation[4] = 5;
	reflectionOperation[5] = 4;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8HBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 3.0e0 * auxMath_._pi / 5.0e0;;
	vector<double> vectorRotations(44);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.000000000000;
	mol0[0].y = 0.000000000000;
	mol0[0].z = 1.000000000000;
	mol0[1].x = 0.000000000000;
	mol0[1].y = 0.000000000000;
	mol0[1].z = -1.000000000000;
	mol0[2].x = 0.000000000000;
	mol0[2].y = -1.000000000000;
	mol0[2].z = 0.000000000000;
	mol0[3].x = 0.866025403784;
	mol0[3].y = -0.500000000000;	
	mol0[3].z = 0.000000000000;
	mol0[4].x = 0.866025403784;
	mol0[4].y = 0.500000000000;
	mol0[4].z = 0.000000000000;
	mol0[5].x = 0.000000000000;
	mol0[5].y = 1.000000000000;
	mol0[5].z = 0.000000000000;
	mol0[6].x = -0.866025403784;
	mol0[6].y = 0.500000000000;
	mol0[6].z = 0.000000000000;
	mol0[7].x = -0.866025403784;
	mol0[7].y = -0.500000000000;	
	mol0[7].z = 0.000000000000;

	//C6-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 3.0e0;
	//C3-1-p
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C2-1
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = auxMath_._pi;
	//C3-1-m
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = -2.0e0 * auxMath_._pi / 3.0e0;
	//C6-1-m
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = -auxMath_._pi / 3.0e0;
	//C2-2 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	vectorRotations[20] = mol0[2].x;
	vectorRotations[21] = mol0[2].y;
	vectorRotations[22] = mol0[2].z;
	vectorRotations[23] = auxMath_._pi;
	//C2-3 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	vectorRotations[24] = mol0[3].x;
	vectorRotations[25] = mol0[3].y;
	vectorRotations[26] = mol0[3].z;
	vectorRotations[27] = auxMath_._pi;
	//C2-4 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	vectorRotations[28] = mol0[4].x;
	vectorRotations[29] = mol0[4].y;
	vectorRotations[30] = mol0[4].z;
	vectorRotations[31] = auxMath_._pi;
	//C2-5 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = auxMath_._pi;
	//C2-6 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[36] = auxReferenceAxis[0];
	vectorRotations[37] = auxReferenceAxis[1];
	vectorRotations[38] = auxReferenceAxis[2];
	vectorRotations[39] = auxMath_._pi;
	//C2-7 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[3].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[3].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[3].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[40] = auxReferenceAxis[0];
	vectorRotations[41] = auxReferenceAxis[1];
	vectorRotations[42] = auxReferenceAxis[2];
	vectorRotations[43] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 1;
	reflectionOperation[1] = 0;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8CU(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(92);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.577350269190;
	mol0[0].y = 0.577350269190;
	mol0[0].z = 0.577350269190;
	mol0[1].x = 0.577350269190;
	mol0[1].y = 0.577350269190; 
	mol0[1].z = -0.577350269190;
	mol0[2].x = 0.577350269190;
	mol0[2].y = -0.577350269190;	
	mol0[2].z = 0.577350269190;
	mol0[3].x = -0.577350269190;
	mol0[3].y = 0.577350269190;	
	mol0[3].z = 0.577350269190;
	mol0[4].x = 0.577350269190;
	mol0[4].y = -0.577350269190;
	mol0[4].z = -0.577350269190;
	mol0[5].x = -0.577350269190;
	mol0[5].y = 0.577350269190;
	mol0[5].z = -0.577350269190;
	mol0[6].x = -0.577350269190;
	mol0[6].y = -0.577350269190;
	mol0[6].z = 0.577350269190;
	mol0[7].x = -0.577350269190;
	mol0[7].y = -0.577350269190;
	mol0[7].z = -0.577350269190;

	//C4-1-p [0-6]
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C2-1 ( C4-2-p C2-2 C4-2-m C4-3-p C2-3 C4-3-m C2-4 C2-7 ) [0-6]
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	//C4-1-m
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = -auxMath_._pi / 2.0e0;

	//C4-2-p
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi / 2.0e0;
	//C2-2 ( C4-1-p C2-1 C4-1-m C4-3-p C2-3 C4-3-m C2-5 C2-8 )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = auxMath_._pi;
	//C4-2-m
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = -auxMath_._pi / 2.0e0;

	//C4-3-p
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = auxMath_._pi / 2.0e0;
	//C2-3 ( C4-1-p C2-1 C4-1-m C4-2-p C2-2 C4-2-m C2-6 C2-9 )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[28] = auxReferenceAxis[0];
	vectorRotations[29] = auxReferenceAxis[1];
	vectorRotations[30] = auxReferenceAxis[2];
	vectorRotations[31] = auxMath_._pi;
	//C4-3-m
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = -auxMath_._pi / 2.0e0;

	//C3-1-p
	vectorRotations[36] = mol0[0].x;
	vectorRotations[37] = mol0[0].y;
	vectorRotations[38] = mol0[0].z;
	vectorRotations[39] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m
	vectorRotations[40] = mol0[0].x;
	vectorRotations[41] = mol0[0].y;
	vectorRotations[42] = mol0[0].z;
	vectorRotations[43] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C3-2-p
	vectorRotations[44] = mol0[1].x;
	vectorRotations[45] = mol0[1].y;
	vectorRotations[46] = mol0[1].z;
	vectorRotations[47] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-2-m
	vectorRotations[48] = mol0[1].x;
	vectorRotations[49] = mol0[1].y;
	vectorRotations[50] = mol0[1].z;
	vectorRotations[51] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C3-3-p
	vectorRotations[52] = mol0[2].x;
	vectorRotations[53] = mol0[2].y;
	vectorRotations[54] = mol0[2].z;
	vectorRotations[55] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-3-m
	vectorRotations[56] = mol0[2].x;
	vectorRotations[57] = mol0[2].y;
	vectorRotations[58] = mol0[2].z;
	vectorRotations[59] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C3-4-p
	vectorRotations[60] = mol0[4].x;
	vectorRotations[61] = mol0[4].y;
	vectorRotations[62] = mol0[4].z;
	vectorRotations[63] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-4-m
	vectorRotations[64] = mol0[4].x;
	vectorRotations[65] = mol0[4].y;
	vectorRotations[66] = mol0[4].z;
	vectorRotations[67] = 4.0e0 * auxMath_._pi / 3.0e0;

	//C2-4 ( C3-3-p C3-3-m C3-4-p C3-4-m C4-1-p C2-1 C4-1-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[68] = auxReferenceAxis[0];
	vectorRotations[69] = auxReferenceAxis[1];
	vectorRotations[70] = auxReferenceAxis[2];
	vectorRotations[71] = auxMath_._pi;
	//C2-5 ( C3-2-p C3-2-m C3-4-p C4-4-m C4-2-p C2-2 C4-2-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[72] = auxReferenceAxis[0];
	vectorRotations[73] = auxReferenceAxis[1];
	vectorRotations[74] = auxReferenceAxis[2];
	vectorRotations[75] = auxMath_._pi;
	//C2-6 ( C3-2-p C3-2-m C3-3-p C3-3-m C4-3-p C2-3 C4-3-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[76] = auxReferenceAxis[0];
	vectorRotations[77] = auxReferenceAxis[1];
	vectorRotations[78] = auxReferenceAxis[2];
	vectorRotations[79] = auxMath_._pi;
	//C2-7 ( C3-1-p C3-1-m C3-2-p C3-2-m C4-1-p C2-1 C4-1-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[2].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[2].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[2].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[80] = auxReferenceAxis[0];
	vectorRotations[81] = auxReferenceAxis[1];
	vectorRotations[82] = auxReferenceAxis[2];
	vectorRotations[83] = auxMath_._pi;
	//C2-8 ( C3-1-p C3-1-m C3-3-p C3-3-m C4-2-p C2-2 C4-2-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[1].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[1].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[1].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[84] = auxReferenceAxis[0];
	vectorRotations[85] = auxReferenceAxis[1];
	vectorRotations[86] = auxReferenceAxis[2];
	vectorRotations[87] = auxMath_._pi;
	//C2-9 ( C3-1-p C3-1-m C3-4-p C3-4-m C4-3-p C2-3 C4-3-m )
	auxReferenceAxis[0] = 0.5e0 *(mol0[1].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[1].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[1].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[88] = auxReferenceAxis[0];
	vectorRotations[89] = auxReferenceAxis[1];
	vectorRotations[90] = auxReferenceAxis[2];
	vectorRotations[91] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[3] = 1;
	reflectionOperation[1] = 3;
	reflectionOperation[4] = 6;
	reflectionOperation[6] = 4;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8ETBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * auxMath_._pi / 180.0e0;;
	vector<double> vectorRotations(20);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.65465368;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.75592894;
	mol0[1].x = -0.65465368;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.75592894;
	mol0[2].x = 0.65465367;
	mol0[2].y = 0.65465367;
	mol0[2].z = -0.37796449;
	mol0[3].x = -0.65465367;
	mol0[3].y = 0.65465367;
	mol0[3].z = -0.37796449;
	mol0[4].x = 0.65465367;
	mol0[4].y = -0.65465367;
	mol0[4].z = -0.37796449;
	mol0[5].x = -0.65465367;
	mol0[5].y = -0.65465367;
	mol0[5].z = -0.37796449;
	mol0[6].x = 1.00000000;
	mol0[6].y = 0.00000000;
	mol0[6].z = 0.00000000;
	mol0[7].x = -1.00000000;
	mol0[7].y = 0.00000000;
	mol0[7].z = 0.00000000;


	//C3-1-p
	vectorRotations[0] = mol0[7].x;
	vectorRotations[1] = mol0[7].y;
	vectorRotations[2] = mol0[7].z;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m
	vectorRotations[4] = mol0[7].x;
	vectorRotations[5] = mol0[7].y;
	vectorRotations[6] = mol0[7].z;
	vectorRotations[7] = -2.0e0 * auxMath_._pi / 3.0e0;

	auxReferenceAxis[0] = mol0[0].x + mol0[1].x;
	auxReferenceAxis[1] = mol0[0].y + mol0[1].y;
	auxReferenceAxis[2] = mol0[0].z + mol0[1].z;
	auxMath_.normalize(auxReferenceAxis);
	//C2-1 ( C3-1-p C3-1-m )
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = auxMath_._pi;

	auxReferenceAxis[0] = mol0[2].x + mol0[3].x;
	auxReferenceAxis[1] = mol0[2].y + mol0[3].y;
	auxReferenceAxis[2] = mol0[2].z + mol0[3].z;
	auxMath_.normalize(auxReferenceAxis);
	//C2-2 ( C3-1-p C3-1-m )
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;

	auxReferenceAxis[0] = mol0[4].x + mol0[5].x;
	auxReferenceAxis[1] = mol0[4].y + mol0[5].y;
	auxReferenceAxis[2] = mol0[4].z + mol0[5].z;
	auxMath_.normalize(auxReferenceAxis);
	//C2-3 ( C3-1-p C3-1-m )
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = auxMath_._pi;


	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[3] = 5;
	reflectionOperation[5] = 3;
	reflectionOperation[2] = 4;
	reflectionOperation[4] = 2;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8HPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * auxMath_._pi / 180.0e0;;
	vector<double> vectorRotations(24);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.62348980;
	mol0[2].y = 0.78183148;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.22252093;
	mol0[3].y = 0.97492791;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.90096887;
	mol0[4].y = 0.43388374;
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.90096887;
	mol0[5].y = -0.43388374;
	mol0[5].z = 0.00000000;
	mol0[6].x = -0.22252093;
	mol0[6].y = -0.97492791;
	mol0[6].z = 0.00000000;
	mol0[7].x = 0.62348980;
	mol0[7].y = -0.78183148;
	mol0[7].z = 0.00000000;

	//C7-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-pp
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-ppp
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-mmm
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-mm
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = 5.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-m
	vectorRotations[20] = mol0[0].x;
	vectorRotations[21] = mol0[0].y;
	vectorRotations[22] = mol0[0].z;
	vectorRotations[23] = 6.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 2;
	reflectionOperation[2] = 1;
	reflectionOperation[3] = 7;
	reflectionOperation[7] = 3;
	reflectionOperation[4] = 6;
	reflectionOperation[6] = 4;

	return vectorRotations;
}



std::vector<double> Geometries::geometry8OP(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * auxMath_._pi / 180.0e0;;
	vector<double> vectorRotations(60);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.70710678;
	mol0[1].y = 0.70710678;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.00000000;
	mol0[2].y = 1.00000000;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.70710678;
	mol0[3].y = 0.70710678;
	mol0[3].z = 0.00000000;
	mol0[4].x = -1.00000000;
	mol0[4].y = 0.00000000;
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.70710678;
	mol0[5].y = -0.70710678;
	mol0[5].z = 0.00000000;
	mol0[6].x = 0.00000000;
	mol0[6].y = -1.00000000;
	mol0[6].z = 0.00000000;
	mol0[7].x = 0.70710678;
	mol0[7].y = -0.70710678;
	mol0[7].z = 0.00000000;

	auxReferenceAxis[0] = 0.0e0;
	auxReferenceAxis[1] = 0.0e0;
	auxReferenceAxis[2] = 1.0e0;

	//C8-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 8.0e0);
	//C4-1-p
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C8-1-ppp
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C2-1 ( C2-2 C2-3 C2-4 C2-5 C2-6 C2-7 C2-8 C2-9 )
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C8-1-mmm
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = 5.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C4-1-m
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = 6.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C8-1-m
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = 7.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);


	//C2-2 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-6 )
	vectorRotations[28] = mol0[0].x;
	vectorRotations[29] = mol0[0].y;
	vectorRotations[30] = mol0[0].z;
	vectorRotations[31] = auxMath_._pi;

	//C2-3 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-7 )
	auxReferenceAxis[0] = mol0[0].x + mol0[1].x;
	auxReferenceAxis[1] = mol0[0].y + mol0[1].y;
	auxReferenceAxis[2] = mol0[0].z + mol0[1].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = auxMath_._pi;

	//C2-4 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-8 )
	vectorRotations[36] = mol0[1].x;
	vectorRotations[37] = mol0[1].y;
	vectorRotations[38] = mol0[1].z;
	vectorRotations[39] = auxMath_._pi;

	//C2-5 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-9 )
	auxReferenceAxis[0] = mol0[1].x + mol0[2].x;
	auxReferenceAxis[1] = mol0[1].y + mol0[2].y;
	auxReferenceAxis[2] = mol0[1].z + mol0[2].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[40] = auxReferenceAxis[0];
	vectorRotations[41] = auxReferenceAxis[1];
	vectorRotations[42] = auxReferenceAxis[2];
	vectorRotations[43] = auxMath_._pi;

	//C2-6 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-2 )
	vectorRotations[44] = mol0[2].x;
	vectorRotations[45] = mol0[2].y;
	vectorRotations[46] = mol0[2].z;
	vectorRotations[47] = auxMath_._pi;

	//C2-7 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-3 )
	auxReferenceAxis[0] = mol0[2].x + mol0[3].x;
	auxReferenceAxis[1] = mol0[2].y + mol0[3].y;
	auxReferenceAxis[2] = mol0[2].z + mol0[3].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[48] = auxReferenceAxis[0];
	vectorRotations[49] = auxReferenceAxis[1];
	vectorRotations[50] = auxReferenceAxis[2];
	vectorRotations[51] = auxMath_._pi;

	//C2-8 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-4 )
	vectorRotations[52] = mol0[3].x;
	vectorRotations[53] = mol0[3].y;
	vectorRotations[54] = mol0[3].z;
	vectorRotations[55] = auxMath_._pi;

	//C2-9 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-5 )
	auxReferenceAxis[0] = mol0[3].x + mol0[4].x;
	auxReferenceAxis[1] = mol0[3].y + mol0[4].y;
	auxReferenceAxis[2] = mol0[3].z + mol0[4].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[56] = auxReferenceAxis[0];
	vectorRotations[57] = auxReferenceAxis[1];
	vectorRotations[58] = auxReferenceAxis[2];
	vectorRotations[59] = auxMath_._pi;


	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;

	return vectorRotations;
}


std::vector<double> Geometries::geometry8JGBF(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 8;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * auxMath_._pi / 180.0e0;;
	vector<double> vectorRotations(12);
	vector<double> auxReferenceAxis(3);
	
	mol0[0].x = 0.50000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.86602540;
	mol0[1].x = -0.50000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.86602540;
	mol0[2].x = 0.70710678;
	mol0[2].y = 0.70710678;
	mol0[2].z = 0.00000000;
	mol0[3].x = 0.70710678;
	mol0[3].y = -0.70710678;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.70710678;
	mol0[4].y = -0.70710678;
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.70710678;
	mol0[5].y = 0.70710678;
	mol0[5].z = 0.00000000;
	mol0[6].x = 0.00000000;
	mol0[6].y = 0.50000000;
	mol0[6].z = -0.86602540;
	mol0[7].x = 0.00000000;
	mol0[7].y = -0.50000000;
	mol0[7].z = -0.86602540;

	//C2-1 ( C2-2 C2-3 )
	auxReferenceAxis[0] = mol0[0].x + mol0[1].x;
	auxReferenceAxis[1] = mol0[0].y + mol0[1].y;
	auxReferenceAxis[2] = mol0[0].z + mol0[1].z;
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;
	//C2-2 ( C2-1 C2-3 )
	vectorRotations[4] = mol0[2].x;
	vectorRotations[5] = mol0[2].y;
	vectorRotations[6] = mol0[2].z;
	vectorRotations[7] = auxMath_._pi;
	//C2-3 ( C2-1 C2-2 )
	vectorRotations[8] = mol0[3].x;
	vectorRotations[9] = mol0[3].y;
	vectorRotations[10] = mol0[3].z;
	vectorRotations[11] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 1;
	reflectionOperation[1] = 0;
	reflectionOperation[4] = 3;
	reflectionOperation[3] = 4;
	reflectionOperation[2] = 5;
	reflectionOperation[5] = 2;

	return vectorRotations;
}



std::vector<double> Geometries::geometry9TCTPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(20);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.000000;
	mol0[0].y = 0.000000;
	mol0[0].z = 1.000000;
	mol0[1].x = -0.23570226;
	mol0[1].y = 0.91287093;
	mol0[1].z = 0.33333333;
	mol0[2].x = -0.94280904;
	mol0[2].y = 0.00000000;
	mol0[2].z = 0.33333333;
	mol0[3].x = 0.23570226;
	mol0[3].y = -0.91287093;
	mol0[3].z = 0.33333333;
	mol0[4].x = 0.94280904;
	mol0[4].y = 0.00000000;
	mol0[4].z = 0.33333333;
	mol0[5].x = 0.53033009;
	mol0[5].y = 0.68465320;
	mol0[5].z = -0.50000000;
	mol0[6].x = -0.53033009;
	mol0[6].y = -0.68465320;
	mol0[6].z = -0.50000000;
	mol0[7].x = -0.58925565;
	mol0[7].y = 0.45643546;
	mol0[7].z = -0.66666667;
	mol0[8].x = 0.58925565;
	mol0[8].y = -0.45643546;
	mol0[8].z = -0.66666667;

	//C2-1 ( C3-1-p C3-1-m )
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi;
	//C2-2 ( C3-1-p C3-1-m )
	vectorRotations[4] = mol0[6].x;
	vectorRotations[5] = mol0[6].y;
	vectorRotations[6] = mol0[6].z;
	vectorRotations[7] = auxMath_._pi;
	//C2-3 ( C3-1-p C3-1-m )
	vectorRotations[8] = mol0[5].x;
	vectorRotations[9] = mol0[5].y;
	vectorRotations[10] = mol0[5].z;
	vectorRotations[11] = auxMath_._pi;
	// produto vetorial
	auxReferenceAxis = auxMath_.vectorProduct(
		mol0[0].x, mol0[0].y, mol0[0].z,
		mol0[6].x, mol0[6].y, mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	//C3-1-p
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 2.0e0 * auxMath_._pi / 3.0e0;
	//C3-1-m
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = 4.0e0 * auxMath_._pi / 3.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[5] = 0;
	reflectionOperation[0] = 5;
	reflectionOperation[8] = 3;
	reflectionOperation[3] = 8;
	reflectionOperation[7] = 2;
	reflectionOperation[2] = 7;

	return vectorRotations;
}


std::vector<double> Geometries::geometry9MFF(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations;

	mol0[0].x = 0.00000000;
	mol0[0].y = 0.98769232;
	mol0[0].z = 0.20656801;
	mol0[1].x = 0.93914925;
	mol0[1].y = 0.30531271;
	mol0[1].z = 0.20656540;
	mol0[2].x = 0.58040926;
	mol0[2].y = -0.79864099;
	mol0[2].z = 0.20655749;
	mol0[3].x = -0.58040926;
	mol0[3].y = -0.79864099;
	mol0[3].z = 0.20655749;
	mol0[4].x = -0.93914925;
	mol0[4].y = 0.30531271;
	mol0[4].z = 0.20656540;
	mol0[5].x = -0.57992608;
	mol0[5].y = -0.33542691;
	mol0[5].z = -0.69361921;
	mol0[6].x = 0.57992608;
	mol0[6].y = -0.33542691;
	mol0[6].z = -0.69361921;
	mol0[7].x = 0.00000000;
	mol0[7].y = 0.66958376;
	mol0[7].z = -0.69426154;
	mol0[8].x = 0.00000000;
	mol0[8].y = 0.00023431;
	mol0[8].z = 1.04868617;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 4;
	reflectionOperation[4] = 1;
	reflectionOperation[2] = 3;
	reflectionOperation[3] = 2;
	reflectionOperation[5] = 6;
	reflectionOperation[6] = 5;

	return vectorRotations;
}



std::vector<double> Geometries::geometry9CSAPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(12);

	mol0[0].x = 0.0000;
	mol0[0].y = 0.0000;
	mol0[0].z = 1.0000;
	mol0[1].x = 0.9322;
	mol0[1].y = 0.0000;
	mol0[1].z = 0.3619;
	mol0[2].x = 0.0000;
	mol0[2].y = 0.9322;
	mol0[2].z = 0.3619;
	mol0[3].x = -0.9322;
	mol0[3].y = 0.0000;
	mol0[3].z = 0.3619;
	mol0[4].x = 0.0000;
	mol0[4].y = -0.9322;
	mol0[4].z = 0.3619;
	mol0[5].x = 0.5606;
	mol0[5].y = 0.5606;
	mol0[5].z = -0.6095;
	mol0[6].x = -0.5606;
	mol0[6].y = 0.5606;
	mol0[6].z = -0.6095;
	mol0[7].x = -0.5606;
	mol0[7].y = -0.5606;
	mol0[7].z = -0.6095;
	mol0[8].x = 0.5606;
	mol0[8].y = -0.5606;
	mol0[8].z = -0.6095;

	//C4-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C2-1
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = auxMath_._pi;
	//C4-1-m
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = -auxMath_._pi / 2.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 4;
	reflectionOperation[4] = 1;
	reflectionOperation[2] = 3;
	reflectionOperation[3] = 2;
	reflectionOperation[5] = 7;
	reflectionOperation[7] = 5;

	return vectorRotations;
}


std::vector<double> Geometries::geometry9CCU(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(12);

	mol0[0].x = 0.64183300;
	mol0[0].y = 0.64183300;
	mol0[0].z = 0.41964300;
	mol0[1].x = 0.64183300;
	mol0[1].y = -0.64183300;
	mol0[1].z = 0.41964300;
	mol0[2].x = -0.64183300;
	mol0[2].y = 0.64183300;
	mol0[2].z = 0.41964300;
	mol0[3].x = -0.64183300;
	mol0[3].y = -0.64183300;
	mol0[3].z = 0.41964300;
	mol0[4].x = 0.53868200;
	mol0[4].y = 0.53868200;
	mol0[4].z = -0.64779900;
	mol0[5].x = 0.53868200;
	mol0[5].y = -0.53868200;
	mol0[5].z = -0.64779900;
	mol0[6].x = -0.53868200;
	mol0[6].y = 0.53868200;
	mol0[6].z = -0.64779900;
	mol0[7].x = -0.53868200;
	mol0[7].y = -0.53868200;
	mol0[7].z = -0.64779900;
	mol0[8].x = 0.00000000;
	mol0[8].y = 0.00000000;
	mol0[8].z = 1.00000000;

	//C4-1-p
	vectorRotations[0] = mol0[8].x;
	vectorRotations[1] = mol0[8].y;
	vectorRotations[2] = mol0[8].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C2-1
	vectorRotations[4] = mol0[8].x;
	vectorRotations[5] = mol0[8].y;
	vectorRotations[6] = mol0[8].z;
	vectorRotations[7] = auxMath_._pi;
	//C4-1-m
	vectorRotations[8] = mol0[8].x;
	vectorRotations[9] = mol0[8].y;
	vectorRotations[10] = mol0[8].z;
	vectorRotations[11] = -auxMath_._pi / 2.0e0;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 2;
	reflectionOperation[2] = 1;
	reflectionOperation[5] = 6;
	reflectionOperation[6] = 5;

	return vectorRotations;
}


std::vector<double> Geometries::geometry9HH(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(4);

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.50000000;
	mol0[1].y = 0.86602540;
	mol0[1].z = 0.00000000;
	mol0[2].x = -0.50000000;
	mol0[2].y = 0.86602540;
	mol0[2].z = 0.00000000;
	mol0[3].x = -1.00000000;
	mol0[3].y = 0.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.50000000;
	mol0[4].y = -0.86602540;
	mol0[4].z = 0.00000000;
	mol0[5].x = 0.50000000;
	mol0[5].y = -0.86602540;
	mol0[5].z = 0.00000000;
	mol0[6].x = 0.00000000;
	mol0[6].y = 0.00000000;
	mol0[6].z = 1.00000000;
	mol0[7].x = 0.50000000;
	mol0[7].y = 0.00000000;
	mol0[7].z = -0.86602540;
	mol0[8].x = -0.50000000;
	mol0[8].y = 0.00000000;
	mol0[8].z = -0.86602540;

	//C2-1
	vectorRotations[0] = mol0[6].x;
	vectorRotations[1] = mol0[6].y;
	vectorRotations[2] = mol0[6].z;
	vectorRotations[3] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[2] = 4;
	reflectionOperation[4] = 2;
	reflectionOperation[1] = 5;
	reflectionOperation[5] = 1;

	return vectorRotations;
}


std::vector<double> Geometries::geometry9OPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(28);

	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.70710678;
	mol0[2].y = 0.70710678;
	mol0[2].z = 0.00000000;
	mol0[3].x = 0.00000000;
	mol0[3].y = 1.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.70710678;
	mol0[4].y = 0.70710678;
	mol0[4].z = 0.00000000;
	mol0[5].x = -1.00000000;
	mol0[5].y = 0.00000000;
	mol0[5].z = 0.00000000;
	mol0[6].x = -0.70710678;
	mol0[6].y = -0.70710678;
	mol0[6].z = 0.00000000;
	mol0[7].x = 0.00000000;
	mol0[7].y = -1.00000000;
	mol0[7].z = 0.00000000;
	mol0[8].x = 0.70710678;
	mol0[8].y = -0.70710678;
	mol0[8].z = 0.00000000;

	//C8-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 8.0e0);
	//C4-1-p
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C8-1-ppp
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C2-1
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C8-1-m
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = -(2.0e0 * auxMath_._pi / 8.0e0);
	//C4-1-m
	vectorRotations[20] = mol0[0].x;
	vectorRotations[21] = mol0[0].y;
	vectorRotations[22] = mol0[0].z;
	vectorRotations[23] = -2.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);
	//C8-1-mmm
	vectorRotations[24] = mol0[0].x;
	vectorRotations[25] = mol0[0].y;
	vectorRotations[26] = mol0[0].z;
	vectorRotations[27] = -3.0e0 * (2.0e0 * auxMath_._pi / 8.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 7;
	reflectionOperation[7] = 1;
	reflectionOperation[2] = 6;
	reflectionOperation[6] = 2;
	reflectionOperation[3] = 5;
	reflectionOperation[5] = 3;

	return vectorRotations;
}


std::vector<double> Geometries::geometry9EP(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(68);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.76604444;
	mol0[1].y = 0.64278761;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.17364818;
	mol0[2].y = 0.98480775;
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.50000000;
	mol0[3].y = 0.86602540;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.93969262;
	mol0[4].y = 0.34202014;
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.93969262;
	mol0[5].y = -0.34202014;
	mol0[5].z = 0.00000000;
	mol0[6].x = -0.50000000;
	mol0[6].y = -0.86602540;
	mol0[6].z = 0.00000000;
	mol0[7].x = 0.17364818;
	mol0[7].y = -0.98480775;
	mol0[7].z = 0.00000000;
	mol0[8].x = 0.76604444;
	mol0[8].y = -0.64278761;
	mol0[8].z = 0.00000000;

	auxReferenceAxis[0] = 0.0e0;
	auxReferenceAxis[1] = 0.0e0;
	auxReferenceAxis[2] = 1.0e0;

	//C9-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (40.e0 * auxMath_._pi / 180.0e0);
	//C9-1-pp
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (40.e0 * auxMath_._pi / 180.0e0);
	//C3-1-p
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * (40.e0 * auxMath_._pi / 180.0e0);
	//C9-1-ppp
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 4.0e0 * (40.e0 * auxMath_._pi / 180.0e0);
	//C9-1-mmm
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = 5.0e0 * (40.e0 * auxMath_._pi / 180.0e0);
	//C3-1-m
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = 6.0e0 * (40.e0 * auxMath_._pi / 180.0e0);
	//C9-1-mm
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = 7.0e0 * (40.e0 * auxMath_._pi / 180.0e0);
	//C9-1-m
	vectorRotations[28] = auxReferenceAxis[0];
	vectorRotations[29] = auxReferenceAxis[1];
	vectorRotations[30] = auxReferenceAxis[2];
	vectorRotations[31] = 8.0e0 * (40.e0 * auxMath_._pi / 180.0e0);

	//C2-1 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[32] = mol0[0].x;
	vectorRotations[33] = mol0[0].y;
	vectorRotations[34] = mol0[0].z;
	vectorRotations[35] = auxMath_._pi;
	//C2-2 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[36] = mol0[1].x;
	vectorRotations[37] = mol0[1].y;
	vectorRotations[38] = mol0[1].z;
	vectorRotations[39] = auxMath_._pi;
	//C2-3 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[40] = mol0[2].x;
	vectorRotations[41] = mol0[2].y;
	vectorRotations[42] = mol0[2].z;
	vectorRotations[43] = auxMath_._pi;
	//C2-4 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[44] = mol0[3].x;
	vectorRotations[45] = mol0[3].y;
	vectorRotations[46] = mol0[3].z;
	vectorRotations[47] = auxMath_._pi;
	//C2-5 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[48] = mol0[4].x;
	vectorRotations[49] = mol0[4].y;
	vectorRotations[50] = mol0[4].z;
	vectorRotations[51] = auxMath_._pi;
	//C2-6 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[52] = mol0[5].x;
	vectorRotations[53] = mol0[5].y;
	vectorRotations[54] = mol0[5].z;
	vectorRotations[55] = auxMath_._pi;
	//C2-7 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[56] = mol0[6].x;
	vectorRotations[57] = mol0[6].y;
	vectorRotations[58] = mol0[6].z;
	vectorRotations[59] = auxMath_._pi;
	//C2-8 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[60] = mol0[7].x;
	vectorRotations[61] = mol0[7].y;
	vectorRotations[62] = mol0[7].z;
	vectorRotations[63] = auxMath_._pi;
	//C2-9 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )
	vectorRotations[64] = mol0[8].x;
	vectorRotations[65] = mol0[8].y;
	vectorRotations[66] = mol0[8].z;
	vectorRotations[67] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;

	return vectorRotations;
}

std::vector<double> Geometries::geometry9HBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(52);

	mol0[0].x = 0.00000000;	
	mol0[0].y = 0.00000000;
	mol0[0].z = -1.00000000;
	mol0[1].x = 1.00000000;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.00000000;
	mol0[2].x = 0.62348980;
	mol0[2].y = 0.78183148;	
	mol0[2].z = 0.00000000;
	mol0[3].x = -0.22252093;
	mol0[3].y = 0.97492791;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.90096887;
	mol0[4].y = 0.43388374;	
	mol0[4].z = 0.00000000;
	mol0[5].x = -0.90096887;
	mol0[5].y = -0.43388374;
	mol0[5].z = 0.00000000;
	mol0[6].x = -0.22252093;
	mol0[6].y = -0.97492791;
	mol0[6].z = 0.00000000;
	mol0[7].x = 0.62348980;
	mol0[7].y = -0.78183148;
	mol0[7].z = 0.00000000;
	mol0[8].x = 0.00000000;
	mol0[8].y = 0.00000000;
	mol0[8].z = 1.00000000;


	//C7-1-p
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-pp
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-ppp
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-mmm
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-mm
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = 5.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);
	//C7-1-m
	vectorRotations[20] = mol0[0].x;
	vectorRotations[21] = mol0[0].y;
	vectorRotations[22] = mol0[0].z;
	vectorRotations[23] = 6.0e0 * (2.0e0 * auxMath_._pi / 7.0e0);


	//C2-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[24] = mol0[1].x;
	vectorRotations[25] = mol0[1].y;
	vectorRotations[26] = mol0[1].z;
	vectorRotations[27] = auxMath_._pi;
	//C2-2 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[28] = mol0[2].x;
	vectorRotations[29] = mol0[2].y;
	vectorRotations[30] = mol0[2].z;
	vectorRotations[31] = auxMath_._pi;
	//C2-3 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[32] = mol0[3].x;
	vectorRotations[33] = mol0[3].y;
	vectorRotations[34] = mol0[3].z;
	vectorRotations[35] = auxMath_._pi;
	//C2-4 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[36] = mol0[4].x;
	vectorRotations[37] = mol0[4].y;
	vectorRotations[38] = mol0[4].z;
	vectorRotations[39] = auxMath_._pi;
	//C2-5 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[40] = mol0[5].x;
	vectorRotations[41] = mol0[5].y;
	vectorRotations[42] = mol0[5].z;
	vectorRotations[43] = auxMath_._pi;
	//C2-6 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[44] = mol0[6].x;
	vectorRotations[45] = mol0[6].y;
	vectorRotations[46] = mol0[6].z;
	vectorRotations[47] = auxMath_._pi;
	//C2-7 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	vectorRotations[48] = mol0[7].x;
	vectorRotations[49] = mol0[7].y;
	vectorRotations[50] = mol0[7].z;
	vectorRotations[51] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 8;
	reflectionOperation[8] = 0;

	return vectorRotations;
}


std::vector<double> Geometries::geometry9JTC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(8);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;	
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.50000000;
	mol0[1].y = 0.86602540;
	mol0[1].z = 0.00000000;
	mol0[2].x = -0.50000000;
	mol0[2].y = 0.86602540;
	mol0[2].z = 0.00000000;
	mol0[3].x = -1.00000000;	
	mol0[3].y = 0.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.50000000;
	mol0[4].y = -0.86602540;
	mol0[4].z = 0.00000000;
	mol0[5].x = 0.50000000;
	mol0[5].y = -0.86602540;
	mol0[5].z = 0.00000000;
	mol0[6].x = 0.50000000;
	mol0[6].y = 0.28867513;
	mol0[6].z = -0.81649658;
	mol0[7].x = -0.50000000;
	mol0[7].y = 0.28867513;
	mol0[7].z = -0.81649658;
	mol0[8].x = 0.00000000;
	mol0[8].y = -0.57735027;
	mol0[8].z = -0.81649658;

	auxReferenceAxis[0] = (mol0[6].x + mol0[7].x + mol0[8].x);
	auxReferenceAxis[1] = (mol0[6].y + mol0[7].y + mol0[8].y);
	auxReferenceAxis[2] = (mol0[6].z + mol0[7].z + mol0[8].z);
	auxMath_.normalize(auxReferenceAxis);

	//C3-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 3.0e0);
	//C3-1-m
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 3.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[4] = 5;
	reflectionOperation[5] = 4;
	reflectionOperation[0] = 3;
	reflectionOperation[3] = 0;
	reflectionOperation[6] = 7;
	reflectionOperation[7] = 6;
	reflectionOperation[1] = 2;
	reflectionOperation[2] = 1;

	return vectorRotations;
}

std::vector<double> Geometries::geometry9JTDIC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 9;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 91.0e0 * (auxMath_._pi / 180.0e0);
	vector<double> vectorRotations(8);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = -0.27639320;
	mol0[0].y = 0.85065081;
	mol0[0].z = -0.44721359;
	mol0[1].x = -0.89442719;
	mol0[1].y = 0.00000000;
	mol0[1].z = -0.44721360;
	mol0[2].x = -0.27639320;
	mol0[2].y = -0.85065081;
	mol0[2].z = -0.44721359;
	mol0[3].x = 0.72360680;
	mol0[3].y = -0.52573111;
	mol0[3].z = -0.44721360;
	mol0[4].x = 0.89442719;
	mol0[4].y = 0.00000000;
	mol0[4].z = 0.44721360;
	mol0[5].x = 0.27639320;
	mol0[5].y = 0.85065081;
	mol0[5].z = 0.44721359;
	mol0[6].x = -0.72360680;
	mol0[6].y = -0.52573111;
	mol0[6].z = 0.44721360;
	mol0[7].x = 0.00000000;
	mol0[7].y = 0.00000000;
	mol0[7].z = -1.00000000;
	mol0[8].x = 0.00000000;
	mol0[8].y = 0.00000000;
	mol0[8].z = 1.00000000;

	auxReferenceAxis[0] = (mol0[1].x + mol0[2].x + mol0[7].x);
	auxReferenceAxis[1] = (mol0[1].y + mol0[2].y + mol0[7].y);
	auxReferenceAxis[2] = (mol0[1].z + mol0[2].z + mol0[7].z);
	auxMath_.normalize(auxReferenceAxis);

	//C3-1-p
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 3.0e0);
	//C3-1-m
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 3.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[1] = 7;
	reflectionOperation[7] = 1;
	reflectionOperation[4] = 8;
	reflectionOperation[8] = 4;
	reflectionOperation[3] = 6;
	reflectionOperation[6] = 3;

	return vectorRotations;
}



std::vector<double> Geometries::geometry10PointSphere(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 10;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(4);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.000000;
	mol0[0].y = 0.000000;
	mol0[0].z = 1.000000;
	mol0[1].x = 0.91458473;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.40439433;
	mol0[2].x = 0.26335401;
	mol0[2].y = 0.87584810;
	mol0[2].z = 0.40439432;
	mol0[3].x = -0.76291954;
	mol0[3].y = 0.50439965;
	mol0[3].z = 0.40439433;
	mol0[4].x = -0.38787671;
	mol0[4].y = -0.82826136;
	mol0[4].z = 0.40439433;
	mol0[5].x = 0.52670802;
	mol0[5].y = -0.82826136;
	mol0[5].z = -0.19121135;
	mol0[6].x = -0.89351054;
	mol0[6].y = -0.25145533;
	mol0[6].z = -0.37203377;
	mol0[7].x = 0.67837321;
	mol0[7].y = 0.50439965;
	mol0[7].z = -0.53421979;
	mol0[8].x = -0.39464230;
	mol0[8].y = 0.69067213;
	mol0[8].z = -0.60599460;
	mol0[9].x = 0.02107419;
	mol0[9].y = -0.25145533;
	mol0[9].z = -0.96763944;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[9] = 6;
	reflectionOperation[6] = 9;
	reflectionOperation[5] = 4;
	reflectionOperation[4] = 5;
	reflectionOperation[7] = 3;
	reflectionOperation[3] = 7;
	reflectionOperation[0] = 1;
	reflectionOperation[1] = 0;

	//c2 - 1
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[7].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[7].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[7].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;

	return vectorRotations;
}

std::vector<double> Geometries::geometry10TD(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 10;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 1.9e0;;
	vector<double> vectorRotations(4);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = -0.50000000;
	mol0[0].y = 0.86600000;
	mol0[0].z = 0.00000000;
	mol0[1].x = 0.50000000;
	mol0[1].y = 0.86600000;
	mol0[1].z = 0.00000000;
	mol0[2].x = -1.00000000;
	mol0[2].y = 0.00000000;
	mol0[2].z = 0.00000000;
	mol0[3].x = 1.00000000;
	mol0[3].y = 0.00000000;
	mol0[3].z = 0.00000000;
	mol0[4].x = -0.50000000;
	mol0[4].y = -0.86600000;
	mol0[4].z = 0.00000000;
	mol0[5].x = 0.50000000;
	mol0[5].y = -0.86600000;
	mol0[5].z = 0.00000000;
	mol0[6].x = -0.50000000;
	mol0[6].y = 0.00000000;
	mol0[6].z = 0.86600000;
	mol0[7].x = 0.50000000;
	mol0[7].y = 0.00000000;
	mol0[7].z = 0.86600000;
	mol0[8].x = 0.00000000;
	mol0[8].y = 0.50000000;
	mol0[8].z = -0.86600000;
	mol0[9].x = 0.00000000;
	mol0[9].y = -0.50000000;
	mol0[9].z = -0.86600000;

	//c2 - 1
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[7].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[7].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[7].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 1;
	reflectionOperation[1] = 0;
	reflectionOperation[2] = 3;
	reflectionOperation[3] = 2;
	reflectionOperation[4] = 5;
	reflectionOperation[5] = 4;
	reflectionOperation[6] = 7;
	reflectionOperation[7] = 6;


	return vectorRotations;
}


std::vector<double> Geometries::geometry10JSPC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 10;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(4);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = -1.72281000;
	mol0[0].y = -0.14410000;
	mol0[0].z = -0.99934000;
	mol0[1].x = -1.72311000;
	mol0[1].y = -0.13180000;
	mol0[1].z = 1.00056000;
	mol0[2].x = -0.88791000;
	mol0[2].y = 1.37940000;
	mol0[2].z = -0.00864000;
	mol0[3].x = 0.04939000;
	mol0[3].y = 0.57650000;
	mol0[3].z = -1.58244000;
	mol0[4].x = -0.11061000;
	mol0[4].y = -1.32760000;
	mol0[4].z = -0.99184000;
	mol0[5].x = -0.11091000;
	mol0[5].y = -1.31520000;
	mol0[5].z = 1.00816000;
	mol0[6].x = 0.04889000;
	mol0[6].y = 0.59600000;
	mol0[6].z = 1.57516000;
	mol0[7].x = 1.10509000;
	mol0[7].y = 1.21240000;
	mol0[7].z = -0.00734000;
	mol0[8].x = 1.67609000;
	mol0[8].y = -0.42900000;
	mol0[8].z = -0.99714000;
	mol0[9].x = 1.67589000;
	mol0[9].y = -0.41660000;
	mol0[9].z = 1.00286000;

	//c2 - 1
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 8;
	reflectionOperation[8] = 0;
	reflectionOperation[1] = 9;
	reflectionOperation[9] = 1;
	reflectionOperation[2] = 7;
	reflectionOperation[7] = 2;

	return vectorRotations;
}

std::vector<double> Geometries::geometry10JBCSAPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 10;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(28);
	vector<double> auxReferenceAxis(3);

	mol0[0].x = 1.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 0.59460000;
	mol0[1].x = 0.70710000;
	mol0[1].y = 0.70710000;
	mol0[1].z = -0.59460000;
	mol0[2].x = 0.00000000;
	mol0[2].y = 1.00000000;
	mol0[2].z = 0.59460000;
	mol0[3].x = -0.70710000;
	mol0[3].y = 0.70710000;
	mol0[3].z = -0.59460000;
	mol0[4].x = -1.00000000;
	mol0[4].y = 0.00000000;
	mol0[4].z = 0.59460000;
	mol0[5].x = -0.70710000;
	mol0[5].y = -0.70710000;
	mol0[5].z = -0.59460000;
	mol0[6].x = 0.00000000;
	mol0[6].y = -1.00000000;
	mol0[6].z = 0.59460000;
	mol0[7].x = 0.70710000;
	mol0[7].y = -0.70710000;
	mol0[7].z = -0.59460000;
	mol0[8].x = 0.00000000;
	mol0[8].y = 0.00000000;
	mol0[8].z = 1.59460000;
	mol0[9].x = 0.00000000;
	mol0[9].y = 0.00000000;
	mol0[9].z = -1.59460000;

	//C4 - 1
	vectorRotations[0] = mol0[9].x;
	vectorRotations[1] = mol0[9].y;
	vectorRotations[2] = mol0[9].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//C4 - 2
	vectorRotations[4] = mol0[9].x;
	vectorRotations[5] = mol0[9].y;
	vectorRotations[6] = mol0[9].z;
	vectorRotations[7] = auxMath_._pi;
	//C4 - 3
	vectorRotations[8] = mol0[9].x;
	vectorRotations[9] = mol0[9].y;
	vectorRotations[10] = mol0[9].z;
	vectorRotations[11] = -auxMath_._pi / 2.0e0;

	//c2 - 1
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;
	//c2 - 2
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = auxMath_._pi;
	//c2 - 3
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = auxMath_._pi;
	//c2 - 4
	auxReferenceAxis.resize(3);
	auxReferenceAxis[0] = 0.5e0 * (mol0[3].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[3].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[3].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[0] = 4;
	reflectionOperation[4] = 0;
	reflectionOperation[1] = 3;
	reflectionOperation[3] = 1;
	reflectionOperation[5] = 7;
	reflectionOperation[7] = 5;

	return vectorRotations;
}

std::vector<double> Geometries::geometry11JCPAPR(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 11;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(16);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.000000;
	mol0[0].y = 0.000000;
	mol0[0].z = 1.000000;
	mol0[1].x = 0.89442719;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.44721360;
	mol0[2].x = 0.27639320;
	mol0[2].y = 0.85065081;
	mol0[2].z = 0.44721360;
	mol0[3].x = -0.72360680;
	mol0[3].y = 0.52573111;
	mol0[3].z = 0.44721360;
	mol0[4].x = -0.72360680;
	mol0[4].y = -0.52573111;
	mol0[4].z = 0.44721360;
	mol0[5].x = 0.27639320;
	mol0[5].y = -0.85065081;
	mol0[5].z = 0.44721360;
	mol0[6].x = 0.72360680;
	mol0[6].y = 0.52573111;
	mol0[6].z = -0.44721360;
	mol0[7].x = -0.27639320;
	mol0[7].y = 0.85065081;
	mol0[7].z = -0.44721360;
	mol0[8].x = -0.89442719;
	mol0[8].y = 0.00000000;
	mol0[8].z = -0.44721360;
	mol0[9].x = -0.27639320;
	mol0[9].y = -0.85065081;
	mol0[9].z = -0.44721360;
	mol0[10].x = 0.00000000;
	mol0[10].y = 0.00000000;
	mol0[10].z = -1.00000000;

	//add rotations
	//c5 - 1
	vectorRotations[0] = mol0[3].x;
	vectorRotations[1] = mol0[3].y;
	vectorRotations[2] = mol0[3].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 5.0e0);
	//c5 - 2
	vectorRotations[4] = mol0[3].x;
	vectorRotations[5] = mol0[3].y;
	vectorRotations[6] = mol0[3].z;
	vectorRotations[7] = 2.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
	//c5 - 3
	vectorRotations[8] = mol0[3].x;
	vectorRotations[9] = mol0[3].y;
	vectorRotations[10] = mol0[3].z;
	vectorRotations[11] = 3.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);
	//c5 - 4
	vectorRotations[12] = mol0[3].x;
	vectorRotations[13] = mol0[3].y;
	vectorRotations[14] = mol0[3].z;
	vectorRotations[15] = 4.0e0 * (2.0e0 * auxMath_._pi / 5.0e0);

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[7] = 2;
	reflectionOperation[2] = 7;
	reflectionOperation[8] = 0;
	reflectionOperation[0] = 8;
	reflectionOperation[10] = 1;
	reflectionOperation[1] = 10;
	reflectionOperation[9] = 5;
	reflectionOperation[5] = 9;


	return vectorRotations;
}

std::vector<double> Geometries::geometry12IC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 12;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
	vector<double> vectorRotations(236);
	vector<double> auxReferenceAxis(3);


	mol0[0].x = 0.00000000;
	mol0[0].y = 0.00000000;
	mol0[0].z = 1.00000000;
	mol0[1].x = 0.89442719;
	mol0[1].y = 0.00000000;
	mol0[1].z = 0.44721360;
	mol0[2].x = 0.27639320;
	mol0[2].y = 0.85065081;
	mol0[2].z = 0.44721360;
	mol0[3].x = -0.72360680;
	mol0[3].y = 0.52573111;
	mol0[3].z = 0.44721360;
	mol0[4].x = -0.72360680;
	mol0[4].y = -0.52573111;
	mol0[4].z = 0.44721360;
	mol0[5].x = 0.27639320;
	mol0[5].y = -0.85065081;
	mol0[5].z = 0.44721360;
	mol0[6].x = 0.72360680;
	mol0[6].y = 0.52573111;
	mol0[6].z = -0.44721360;
	mol0[7].x = -0.27639320;
	mol0[7].y = 0.85065081;
	mol0[7].z = -0.44721360;
	mol0[8].x = -0.89442719;
	mol0[8].y = 0.00000000;
	mol0[8].z = -0.44721360;
	mol0[9].x = -0.27639320;
	mol0[9].y = -0.85065081;
	mol0[9].z = -0.44721360;
	mol0[10].x = 0.72360680;
	mol0[10].y = -0.52573111;
	mol0[10].z = -0.44721360;
	mol0[11].x = 0.00000000;
	mol0[11].y = 0.00000000;
	mol0[11].z = -1.00000000;

	// 4 c5 para cada um dos 6 pontos principais (12 c5 e 12 c52)
	// 2 c3 pra cada face unica (20 c3)
	// c2 para cada aresta (10 no topo e 5 ao redor = 15)
	//add rotations
	// 0 por cima
	// c5 - 0
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = (2.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = (4.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = (6.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = (8.0e0 * auxMath_._pi / 5.0e0);

	// c5 - 1
	vectorRotations[16] = mol0[1].x;
	vectorRotations[17] = mol0[1].y;
	vectorRotations[18] = mol0[1].z;
	vectorRotations[19] = (2.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[20] = mol0[1].x;
	vectorRotations[21] = mol0[1].y;
	vectorRotations[22] = mol0[1].z;
	vectorRotations[23] = (4.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[24] = mol0[1].x;
	vectorRotations[25] = mol0[1].y;
	vectorRotations[26] = mol0[1].z;
	vectorRotations[27] = (6.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[28] = mol0[1].x;
	vectorRotations[29] = mol0[1].y;
	vectorRotations[30] = mol0[1].z;
	vectorRotations[31] = (8.0e0 * auxMath_._pi / 5.0e0);

	// c5 - 2
	vectorRotations[32] = mol0[2].x;
	vectorRotations[33] = mol0[2].y;
	vectorRotations[34] = mol0[2].z;
	vectorRotations[35] = (2.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[36] = mol0[2].x;
	vectorRotations[37] = mol0[2].y;
	vectorRotations[38] = mol0[2].z;
	vectorRotations[39] = (4.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[40] = mol0[2].x;
	vectorRotations[41] = mol0[2].y;
	vectorRotations[42] = mol0[2].z;
	vectorRotations[43] = (6.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[44] = mol0[2].x;
	vectorRotations[45] = mol0[2].y;
	vectorRotations[46] = mol0[2].z;
	vectorRotations[47] = (8.0e0 * auxMath_._pi / 5.0e0);

	// c5 - 3
	vectorRotations[48] = mol0[3].x;
	vectorRotations[49] = mol0[3].y;
	vectorRotations[50] = mol0[3].z;
	vectorRotations[51] = (2.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[52] = mol0[3].x;
	vectorRotations[53] = mol0[3].y;
	vectorRotations[54] = mol0[3].z;
	vectorRotations[55] = (4.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[56] = mol0[3].x;
	vectorRotations[57] = mol0[3].y;
	vectorRotations[58] = mol0[3].z;
	vectorRotations[59] = (6.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[60] = mol0[3].x;
	vectorRotations[61] = mol0[3].y;
	vectorRotations[62] = mol0[3].z;
	vectorRotations[63] = (8.0e0 * auxMath_._pi / 5.0e0);

	// c5 - 4
	vectorRotations[64] = mol0[4].x;
	vectorRotations[65] = mol0[4].y;
	vectorRotations[66] = mol0[4].z;
	vectorRotations[67] = (2.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[68] = mol0[4].x;
	vectorRotations[69] = mol0[4].y;
	vectorRotations[70] = mol0[4].z;
	vectorRotations[71] = (4.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[72] = mol0[4].x;
	vectorRotations[73] = mol0[4].y;
	vectorRotations[74] = mol0[4].z;
	vectorRotations[75] = (6.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[76] = mol0[4].x;
	vectorRotations[77] = mol0[4].y;
	vectorRotations[78] = mol0[4].z;
	vectorRotations[79] = (8.0e0 * auxMath_._pi / 5.0e0);

	// c5 - 5
	vectorRotations[80] = mol0[5].x;
	vectorRotations[81] = mol0[5].y;
	vectorRotations[82] = mol0[5].z;
	vectorRotations[83] = (2.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[84] = mol0[5].x;
	vectorRotations[85] = mol0[5].y;
	vectorRotations[86] = mol0[5].z;
	vectorRotations[87] = (4.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[88] = mol0[5].x;
	vectorRotations[89] = mol0[5].y;
	vectorRotations[90] = mol0[5].z;
	vectorRotations[91] = (6.0e0 * auxMath_._pi / 5.0e0);
	vectorRotations[92] = mol0[5].x;
	vectorRotations[93] = mol0[5].y;
	vectorRotations[94] = mol0[5].z;
	vectorRotations[95] = (8.0e0 * auxMath_._pi / 5.0e0);

	auxReferenceAxis.resize(3);
	// c3 pra cada face
	// 0-1-2
	auxReferenceAxis[0] = (mol0[0].x + mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[96] = auxReferenceAxis[0];
	vectorRotations[97] = auxReferenceAxis[1];
	vectorRotations[98] = auxReferenceAxis[2];
	vectorRotations[99] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[100] = auxReferenceAxis[0];
	vectorRotations[101] = auxReferenceAxis[1];
	vectorRotations[102] = auxReferenceAxis[2];
	vectorRotations[103] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 0-2-3
	auxReferenceAxis[0] = (mol0[0].x + mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[104] = auxReferenceAxis[0];
	vectorRotations[105] = auxReferenceAxis[1];
	vectorRotations[106] = auxReferenceAxis[2];
	vectorRotations[107] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[108] = auxReferenceAxis[0];
	vectorRotations[109] = auxReferenceAxis[1];
	vectorRotations[110] = auxReferenceAxis[2];
	vectorRotations[111] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 0-3-4
	auxReferenceAxis[0] = (mol0[0].x + mol0[3].x + mol0[4].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[3].y + mol0[4].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[3].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[112] = auxReferenceAxis[0];
	vectorRotations[113] = auxReferenceAxis[1];
	vectorRotations[114] = auxReferenceAxis[2];
	vectorRotations[115] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[116] = auxReferenceAxis[0];
	vectorRotations[117] = auxReferenceAxis[1];
	vectorRotations[118] = auxReferenceAxis[2];
	vectorRotations[119] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 0-4-5
	auxReferenceAxis[0] = (mol0[0].x + mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[120] = auxReferenceAxis[0];
	vectorRotations[121] = auxReferenceAxis[1];
	vectorRotations[122] = auxReferenceAxis[2];
	vectorRotations[123] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[124] = auxReferenceAxis[0];
	vectorRotations[125] = auxReferenceAxis[1];
	vectorRotations[126] = auxReferenceAxis[2];
	vectorRotations[127] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 0-5-1
	auxReferenceAxis[0] = (mol0[0].x + mol0[5].x + mol0[1].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[5].y + mol0[1].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[5].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[128] = auxReferenceAxis[0];
	vectorRotations[129] = auxReferenceAxis[1];
	vectorRotations[130] = auxReferenceAxis[2];
	vectorRotations[131] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[132] = auxReferenceAxis[0];
	vectorRotations[133] = auxReferenceAxis[1];
	vectorRotations[134] = auxReferenceAxis[2];
	vectorRotations[135] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 3-8-4
	auxReferenceAxis[0] = (mol0[3].x + mol0[8].x + mol0[4].x);
	auxReferenceAxis[1] = (mol0[3].y + mol0[8].y + mol0[4].y);
	auxReferenceAxis[2] = (mol0[3].z + mol0[8].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[136] = auxReferenceAxis[0];
	vectorRotations[137] = auxReferenceAxis[1];
	vectorRotations[138] = auxReferenceAxis[2];
	vectorRotations[139] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[140] = auxReferenceAxis[0];
	vectorRotations[141] = auxReferenceAxis[1];
	vectorRotations[142] = auxReferenceAxis[2];
	vectorRotations[143] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 8-4-9
	auxReferenceAxis[0] = (mol0[8].x + mol0[4].x + mol0[9].x);
	auxReferenceAxis[1] = (mol0[8].y + mol0[4].y + mol0[9].y);
	auxReferenceAxis[2] = (mol0[8].z + mol0[4].z + mol0[9].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[144] = auxReferenceAxis[0];
	vectorRotations[145] = auxReferenceAxis[1];
	vectorRotations[146] = auxReferenceAxis[2];
	vectorRotations[147] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[148] = auxReferenceAxis[0];
	vectorRotations[149] = auxReferenceAxis[1];
	vectorRotations[150] = auxReferenceAxis[2];
	vectorRotations[151] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 4-9-5
	auxReferenceAxis[0] = (mol0[4].x + mol0[9].x + mol0[5].x);
	auxReferenceAxis[1] = (mol0[4].y + mol0[9].y + mol0[5].y);
	auxReferenceAxis[2] = (mol0[4].z + mol0[9].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[152] = auxReferenceAxis[0];
	vectorRotations[153] = auxReferenceAxis[1];
	vectorRotations[154] = auxReferenceAxis[2];
	vectorRotations[155] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[156] = auxReferenceAxis[0];
	vectorRotations[157] = auxReferenceAxis[1];
	vectorRotations[158] = auxReferenceAxis[2];
	vectorRotations[159] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 9-5-10
	auxReferenceAxis[0] = (mol0[9].x + mol0[5].x + mol0[10].x);
	auxReferenceAxis[1] = (mol0[9].y + mol0[5].y + mol0[10].y);
	auxReferenceAxis[2] = (mol0[9].z + mol0[5].z + mol0[10].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[160] = auxReferenceAxis[0];
	vectorRotations[161] = auxReferenceAxis[1];
	vectorRotations[162] = auxReferenceAxis[2];
	vectorRotations[163] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[164] = auxReferenceAxis[0];
	vectorRotations[165] = auxReferenceAxis[1];
	vectorRotations[166] = auxReferenceAxis[2];
	vectorRotations[167] = 4.0e0 * auxMath_._pi / 3.0e0;

	// 5-10-1
	auxReferenceAxis[0] = (mol0[10].x + mol0[5].x + mol0[1].x);
	auxReferenceAxis[1] = (mol0[10].y + mol0[5].y + mol0[1].y);
	auxReferenceAxis[2] = (mol0[10].z + mol0[5].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[168] = auxReferenceAxis[0];
	vectorRotations[169] = auxReferenceAxis[1];
	vectorRotations[170] = auxReferenceAxis[2];
	vectorRotations[171] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[172] = auxReferenceAxis[0];
	vectorRotations[173] = auxReferenceAxis[1];
	vectorRotations[174] = auxReferenceAxis[2];
	vectorRotations[175] = 4.0e0 * auxMath_._pi / 3.0e0;

	// c2 nas arestas
	// 0-1
	auxReferenceAxis[0] = (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[176] = auxReferenceAxis[0];
	vectorRotations[177] = auxReferenceAxis[1];
	vectorRotations[178] = auxReferenceAxis[2];
	vectorRotations[179] = auxMath_._pi;

	// 0-2
	auxReferenceAxis[0] = (mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[180] = auxReferenceAxis[0];
	vectorRotations[181] = auxReferenceAxis[1];
	vectorRotations[182] = auxReferenceAxis[2];
	vectorRotations[183] = auxMath_._pi;

	// 0-3
	auxReferenceAxis[0] = (mol0[0].x + mol0[3].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[3].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[184] = auxReferenceAxis[0];
	vectorRotations[185] = auxReferenceAxis[1];
	vectorRotations[186] = auxReferenceAxis[2];
	vectorRotations[187] = auxMath_._pi;

	// 0-4
	auxReferenceAxis[0] = (mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[188] = auxReferenceAxis[0];
	vectorRotations[189] = auxReferenceAxis[1];
	vectorRotations[190] = auxReferenceAxis[2];
	vectorRotations[191] = auxMath_._pi;

	// 0-5
	auxReferenceAxis[0] = (mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[192] = auxReferenceAxis[0];
	vectorRotations[193] = auxReferenceAxis[1];
	vectorRotations[194] = auxReferenceAxis[2];
	vectorRotations[195] = auxMath_._pi;

	// 1-2
	auxReferenceAxis[0] = (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[196] = auxReferenceAxis[0];
	vectorRotations[197] = auxReferenceAxis[1];
	vectorRotations[198] = auxReferenceAxis[2];
	vectorRotations[199] = auxMath_._pi;

	// 2-3
	auxReferenceAxis[0] = (mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = (mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = (mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[200] = auxReferenceAxis[0];
	vectorRotations[201] = auxReferenceAxis[1];
	vectorRotations[202] = auxReferenceAxis[2];
	vectorRotations[203] = auxMath_._pi;

	// 3-4
	auxReferenceAxis[0] = (mol0[3].x + mol0[4].x);
	auxReferenceAxis[1] = (mol0[3].y + mol0[4].y);
	auxReferenceAxis[2] = (mol0[3].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[204] = auxReferenceAxis[0];
	vectorRotations[205] = auxReferenceAxis[1];
	vectorRotations[206] = auxReferenceAxis[2];
	vectorRotations[207] = auxMath_._pi;

	// 4-5
	auxReferenceAxis[0] = (mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = (mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = (mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[208] = auxReferenceAxis[0];
	vectorRotations[209] = auxReferenceAxis[1];
	vectorRotations[210] = auxReferenceAxis[2];
	vectorRotations[211] = auxMath_._pi;

	// 5-1
	auxReferenceAxis[0] = (mol0[5].x + mol0[1].x);
	auxReferenceAxis[1] = (mol0[5].y + mol0[1].y);
	auxReferenceAxis[2] = (mol0[5].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[212] = auxReferenceAxis[0];
	vectorRotations[213] = auxReferenceAxis[1];
	vectorRotations[214] = auxReferenceAxis[2];
	vectorRotations[215] = auxMath_._pi;

	// 1-6
	auxReferenceAxis[0] = (mol0[1].x + mol0[6].x);
	auxReferenceAxis[1] = (mol0[1].y + mol0[6].y);
	auxReferenceAxis[2] = (mol0[1].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[216] = auxReferenceAxis[0];
	vectorRotations[217] = auxReferenceAxis[1];
	vectorRotations[218] = auxReferenceAxis[2];
	vectorRotations[219] = auxMath_._pi;

	// 6-2
	auxReferenceAxis[0] = (mol0[6].x + mol0[2].x);
	auxReferenceAxis[1] = (mol0[6].y + mol0[2].y);
	auxReferenceAxis[2] = (mol0[6].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[220] = auxReferenceAxis[0];
	vectorRotations[221] = auxReferenceAxis[1];
	vectorRotations[222] = auxReferenceAxis[2];
	vectorRotations[223] = auxMath_._pi;

	// 2-7
	auxReferenceAxis[0] = (mol0[2].x + mol0[7].x);
	auxReferenceAxis[1] = (mol0[2].y + mol0[7].y);
	auxReferenceAxis[2] = (mol0[2].z + mol0[7].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[224] = auxReferenceAxis[0];
	vectorRotations[225] = auxReferenceAxis[1];
	vectorRotations[226] = auxReferenceAxis[2];
	vectorRotations[227] = auxMath_._pi;

	// 7-3
	auxReferenceAxis[0] = (mol0[7].x + mol0[3].x);
	auxReferenceAxis[1] = (mol0[7].y + mol0[3].y);
	auxReferenceAxis[2] = (mol0[7].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[228] = auxReferenceAxis[0];
	vectorRotations[229] = auxReferenceAxis[1];
	vectorRotations[230] = auxReferenceAxis[2];
	vectorRotations[231] = auxMath_._pi;

	// 3-8
	auxReferenceAxis[0] = (mol0[3].x + mol0[8].x);
	auxReferenceAxis[1] = (mol0[3].y + mol0[8].y);
	auxReferenceAxis[2] = (mol0[3].z + mol0[8].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[232] = auxReferenceAxis[0];
	vectorRotations[233] = auxReferenceAxis[1];
	vectorRotations[234] = auxReferenceAxis[2];
	vectorRotations[235] = auxMath_._pi;

	for (size_t i = 0; i < reflectionOperation.size(); i++)
		reflectionOperation[i] = i;
	reflectionOperation[4] = 1;
	reflectionOperation[1] = 4;
	reflectionOperation[3] = 2;
	reflectionOperation[2] = 3;
	reflectionOperation[9] = 10;
	reflectionOperation[10] = 9;
	reflectionOperation[8] = 6;
	reflectionOperation[6] = 8;

	return vectorRotations;
}










void Geometries::geometry4TetrahedronotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 4;
	allReflections.resize(12);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;

	//P-2
	allReflections[1][0] = 2;
	allReflections[1][2] = 0;

	//P-3
	allReflections[2][0] = 3;
	allReflections[2][3] = 0;

	//P-4
	allReflections[3][1] = 2;
	allReflections[3][2] = 1;

	//P-5
	allReflections[4][1] = 3;
	allReflections[4][3] = 1;

	//P-6
	allReflections[5][2] = 3;
	allReflections[5][3] = 2;

	//S4-1-p
	allReflections[6][0] = 3;
	allReflections[6][3] = 1;
	allReflections[6][1] = 2;
	allReflections[6][2] = 0;

	//S4-1-m
	allReflections[7][0] = 2;
	allReflections[7][2] = 1;
	allReflections[7][1] = 3;
	allReflections[7][3] = 0;

	//S4-2-p
	allReflections[8][0] = 1;
	allReflections[8][1] = 2;
	allReflections[8][2] = 3;
	allReflections[8][3] = 0;

	//S4-2-m
	allReflections[9][0] = 3;
	allReflections[9][3] = 2;
	allReflections[9][2] = 1;
	allReflections[9][1] = 0;

	//S4-3-p
	allReflections[10][0] = 2;
	allReflections[10][2] = 3;
	allReflections[10][3] = 1;
	allReflections[10][1] = 0;

	//S4-3-m
	allReflections[11][0] = 1;
	allReflections[11][1] = 3;
	allReflections[11][3] = 2;
	allReflections[11][2] = 0;

}

string Geometries::geometry4TetrahedronSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
		case 1:
			return "C3-1-m";
		case 2:
			return "C3-2-p";
		case 3:
			return "C3-2-m";
		case 4:
			return "C3-3-p";
		case 5:
			return "C3-3-m";
		case 6:
			return "C3-4-p";
		case 7:
			return "C3-4-m";
		case 8:
			return "C2-1 ( C2-2 C2-3 )";
		case 9:
			return "C2-2 ( C2-1 C2-3 )";
		case 10:
			return "C2-3 ( C2-1 C2-2 )";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "S4-1-p";
		case 7:
			return "S4-1-m";
		case 8:
			return "S4-2-p";
		case 9:
			return "S4-2-m";
		case 10:
			return "S4-3-p";
		case 11:
			return "S4-3-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}


void Geometries::geometry4SquareotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 4;
	allReflections.resize(8);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1 ( C4-1-p C2-1 C4-1-m )

	//P-2 ( C2-4 )
	allReflections[1][0] = 2;
	allReflections[1][2] = 0;

	//P-3 ( C2-2 )
	allReflections[2][1] = 3;
	allReflections[2][3] = 1;

	//P-4 ( C2-3 )
	allReflections[3][0] = 3;
	allReflections[3][3] = 0;
	allReflections[3][1] = 2;
	allReflections[3][2] = 1;

	//P-5 ( C2-5 )
	allReflections[4][0] = 1;
	allReflections[4][1] = 0;
	allReflections[4][2] = 3;
	allReflections[4][3] = 2;

	//Inv 
	allReflections[5][0] = 2;
	allReflections[5][2] = 0;
	allReflections[5][1] = 3;
	allReflections[5][3] = 1;

	//S4-1-p
	allReflections[6][0] = 1;
	allReflections[6][1] = 2;
	allReflections[6][2] = 3;
	allReflections[6][3] = 0;

	//S4-1-m
	allReflections[7][0] = 3;
	allReflections[7][3] = 2;
	allReflections[7][2] = 1;
	allReflections[7][1] = 0;
}

string Geometries::geometry4SquareSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-2 ( C4-1-p C2-1 C4-1-m C2-4 )";
		case 1:
			return "C2-3 ( C4-1-p C2-1 C4-1-m C2-5 )";
		case 2:
			return "C2-4 ( C4-1-p C2-1 C4-1-m C2-2 )";
		case 3:
			return "C2-5 ( C4-1-p C2-1 C4-1-m C2-3 )";
		case 4:
			return "C4-1-p";
		case 5:
			return "C2-1 ( C2-2 C2-3 C2-4 C2-5 )";
		case 6:
			return "C4-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C4-1-p C2-1 C4-1-m )";
		case 1:
			return "P-2 ( C2-4 )";
		case 2:
			return "P-3 ( C2-2 )";
		case 3:
			return "P-4 ( C2-3 )";
		case 4:
			return "P-5 ( C2-5 )";
		case 5:
			return "Inv";
		case 6:
			return "S4-1-p";
		case 7:
			return "S4-1-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}


void Geometries::geometry4SSotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 4;
	allReflections.resize(2);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;

	//P-2
	allReflections[1][0] = 3;
	allReflections[1][3] = 0;
}

string Geometries::geometry4SSSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}


void Geometries::geometry4vTBPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 4;
	allReflections.resize(3);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1 
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;

	//P-2
	allReflections[1][2] = 3;
	allReflections[1][3] = 2;

	//P-3
	allReflections[2][1] = 3;
	allReflections[2][3] = 1;
}

string Geometries::geometry4vTBPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
			break;
		case 1:
			return "C3-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}



void Geometries::geometry5SPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 5;
	allReflections.resize(4);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][3] = 4;
	allReflections[0][4] = 3;

	//P-2
	allReflections[1][1] = 3;
	allReflections[1][3] = 1;
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;

	//P-3
	allReflections[2][1] = 2;
	allReflections[2][2] = 1;

	//P-4
	allReflections[3][2] = 3;
	allReflections[3][3] = 2;
	allReflections[3][1] = 4;
	allReflections[3][4] = 1;
}

string Geometries::geometry5SPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C4-1-p";
			break;

		case 1:
			return "C2-1 ";
			break;

		case 2:
			return "C4-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}

void Geometries::geometry5TBPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 5;
	allReflections.resize(6);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1 ( C3-1-p C3-1-m )
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;

	//P-2
	allReflections[1][0] = 3;
	allReflections[1][3] = 0;

	//P-3
	allReflections[2][4] = 3;
	allReflections[2][3] = 4;

	//P-4
	allReflections[3][0] = 4;
	allReflections[3][4] = 0;

	//S3-1-p
	allReflections[4][1] = 2;
	allReflections[4][2] = 1;
	allReflections[4][0] = 3;
	allReflections[4][3] = 4;
	allReflections[4][4] = 0;

	//S3-1-m
	allReflections[5][1] = 2;
	allReflections[5][2] = 1;
	allReflections[5][0] = 4;
	allReflections[5][4] = 3;
	allReflections[5][3] = 0;

}

string Geometries::geometry5TBPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
			break;

		case 1:
			return "C3-1-m";
			break;

		case 2:
			return "C2-1 ( C3-1-p C3-1-m )";
			break;

		case 3:
			return "C2-2 ( C3-1-p C3-1-m )";
			break;

		case 4:
			return "C2-3 ( C3-1-p C3-1-m )";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C3-1-p C3-1-m )";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "S3-1-p";
		case 5:
			return "S3-1-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}


void Geometries::geometry5PPotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 5;
	allReflections.resize(10);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) 

	//P-2
	allReflections[1][2] = 3;
	allReflections[1][3] = 2;
	allReflections[1][1] = 4;
	allReflections[1][4] = 1;

	//P-3
	allReflections[2][4] = 3;
	allReflections[2][3] = 4;
	allReflections[2][0] = 2;
	allReflections[2][2] = 0;

	//P-4
	allReflections[3][0] = 4;
	allReflections[3][4] = 0;
	allReflections[3][1] = 3;
	allReflections[3][3] = 1;

	//P-5
	allReflections[4][0] = 1;
	allReflections[4][1] = 0;
	allReflections[4][2] = 4;
	allReflections[4][4] = 2;

	//P-6
	allReflections[5][1] = 2;
	allReflections[5][2] = 1;
	allReflections[5][0] = 3;
	allReflections[5][3] = 0;

	//S5-1-p
	allReflections[6][0] = 1;
	allReflections[6][1] = 2;
	allReflections[6][2] = 3;
	allReflections[6][3] = 4;
	allReflections[6][4] = 0;

	//S5-1-pp
	allReflections[7][0] = 2;
	allReflections[7][2] = 4;
	allReflections[7][4] = 1;
	allReflections[7][1] = 3;
	allReflections[7][3] = 0;

	//S5-1-mm
	allReflections[8][0] = 3;
	allReflections[8][3] = 1;
	allReflections[8][1] = 4;
	allReflections[8][4] = 2;
	allReflections[8][2] = 0;

	//S5-1-m
	allReflections[9][0] = 4;
	allReflections[9][4] = 3;
	allReflections[9][3] = 2;
	allReflections[9][2] = 1;
	allReflections[9][1] = 0;

}

string Geometries::geometry5PPSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C5-1-p";
			break;
		case 1:
			return "C5-1-pp";
			break;
		case 2:
			return "C5-1-mm";
			break;
		case 3:
			return "C5-1-m";
			break;
		case 4:
			return "C2-1 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) ";
			break;
		case 5:
			return "C2-2 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) ";
			break;
		case 6:
			return "C2-3 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) ";
			break;
		case 7:
			return "C2-4 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) ";
			break;
		case 8:
			return "C2-5 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) ";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C5-1-p C5-1-pp C5-1-mm C5-1-m ) ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "S5-1-p";
		case 7:
			return "S5-1-pp";
		case 8:
			return "S5-1-mm";
		case 9:
			return "S5-1-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}



void Geometries::geometry6TPRotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 6;
	allReflections.resize(6);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
//	reflectionOperation[0] = 4;
//	reflectionOperation[4] = 0;
//	reflectionOperation[1] = 5;
//	reflectionOperation[5] = 1;

	//P-1 ( C3-1-p C3-1-m )
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;
	allReflections[0][2] = 3;
	allReflections[0][3] = 2;
	allReflections[0][4] = 5;
	allReflections[0][5] = 4;

	//P-2
	allReflections[1][0] = 4;
	allReflections[1][4] = 0;
	allReflections[1][1] = 5;
	allReflections[1][5] = 1;

	//P-3
	allReflections[2][2] = 4;
	allReflections[2][4] = 2;
	allReflections[2][3] = 5;
	allReflections[2][5] = 3;

	//P-4
	allReflections[3][0] = 2;
	allReflections[3][2] = 0;
	allReflections[3][1] = 3;
	allReflections[3][3] = 1;

	//S3-1-p
	allReflections[4][0] = 3;
	allReflections[4][3] = 4;
	allReflections[4][4] = 1;
	allReflections[4][1] = 2;
	allReflections[4][2] = 5;
	allReflections[4][5] = 0;

	//S3-1-m
	allReflections[5][0] = 5;
	allReflections[5][5] = 2;
	allReflections[5][2] = 1;
	allReflections[5][1] = 4;
	allReflections[5][4] = 3;
	allReflections[5][3] = 0;

}

string Geometries::geometry6TPRSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
			break;

		case 1:
			return "C3-1-m";
			break;

		case 2:
			return "C2-1 ( C3-1-p C3-1-m )";
			break;

		case 3:
			return "C2-2 ( C3-1-p C3-1-m )";
			break;

		case 4:
			return "C2-3 ( C3-1-p C3-1-m )";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C3-1-p C3-1-m )";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "S3-1-p";
		case 5:
			return "S3-1-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}


void Geometries::geometry6HPotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 6;
	allReflections.resize(12);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	 
	//P-1 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m )

	//P-2 ( C2-5 )
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;
	allReflections[1][1] = 5;
	allReflections[1][5] = 1;

	//P-3 ( C2-6 )
	allReflections[2][2] = 3;
	allReflections[2][3] = 2;
	allReflections[2][1] = 4;
	allReflections[2][4] = 1;
	allReflections[2][0] = 5;
	allReflections[2][5] = 0;

	//P-4 ( C2-7 )
	allReflections[3][1] = 3;
	allReflections[3][3] = 1;
	allReflections[3][0] = 4;
	allReflections[3][4] = 0;

	//P-5 ( C2-2 )
	allReflections[4][1] = 2;
	allReflections[4][2] = 1;
	allReflections[4][0] = 3;
	allReflections[4][3] = 0;
	allReflections[4][4] = 5;
	allReflections[4][5] = 4;

	//P-6 ( C2-3 )
	allReflections[5][0] = 2;
	allReflections[5][2] = 0;
	allReflections[5][3] = 5;
	allReflections[5][5] = 3;

	//P-7 ( C2-4 )
	allReflections[6][0] = 1;
	allReflections[6][1] = 0;
	allReflections[6][2] = 5;
	allReflections[6][5] = 2;
	allReflections[6][4] = 3;
	allReflections[6][3] = 4;

	//S6-1-p
	allReflections[7][0] = 1;
	allReflections[7][1] = 2;
	allReflections[7][2] = 3;
	allReflections[7][3] = 4;
	allReflections[7][4] = 5;
	allReflections[7][5] = 0;

	//S3-1-p
	allReflections[8][0] = 2;
	allReflections[8][2] = 4;
	allReflections[8][4] = 0;
	allReflections[8][1] = 3;
	allReflections[8][3] = 5;
	allReflections[8][5] = 1;

	//S3-1-m
	allReflections[9][0] = 4;
	allReflections[9][4] = 2;
	allReflections[9][2] = 0;
	allReflections[9][1] = 5;
	allReflections[9][5] = 3;
	allReflections[9][3] = 1;

	//S6-1-m
	allReflections[10][0] = 5;
	allReflections[10][5] = 4;
	allReflections[10][4] = 3;
	allReflections[10][3] = 2;
	allReflections[10][2] = 1;
	allReflections[10][1] = 0;

	//Inv
	allReflections[11][0] = 3;
	allReflections[11][3] = 0;
	allReflections[11][1] = 4;
	allReflections[11][4] = 1;
	allReflections[11][2] = 5;
	allReflections[11][5] = 2;

}

string Geometries::geometry6HPSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C6-1-p";
			break;
		case 1:
			return "C3-1-p";
			break;
		case 2:
			return "C2-1 ( C2-2 C2-3 C2-4 C2-5 C2-6 C2-7 ) ";
			break;
		case 3:
			return "C3-1-m";
			break;
		case 4:
			return "C6-1-m";
			break;
		case 5:
			return "C2-2 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-5 ) ";
			break;
		case 6:
			return "C2-3 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-6 ) ";
			break;
		case 7:
			return "C2-4 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-7 ) ";
			break;
		case 8:
			return "C2-5 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-2 ) ";
			break;
		case 9:
			return "C2-6 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-3 ) ";
			break;
		case 10:
			return "C2-7 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m C2-4 ) ";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C6-1-p C3-1-p C2-1 C3-1-m C6-1-m ) ";
		case 1:
			return "P-2 ( C2-5 ) ";
		case 2:
			return "P-3 ( C2-6 ) ";
		case 3:
			return "P-4 ( C2-7 ) ";
		case 4:
			return "P-5 ( C2-2 ) ";
		case 5:
			return "P-6 ( C2-3 ) ";
		case 6:
			return "P-7 ( C2-4 ) ";
		case 7:
			return "S6-1-p";
		case 8:
			return "S3-1-p";
		case 9:
			return "S3-1-m";
		case 10:
			return "S6-1-m";
		case 11:
			return "Inv";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}



void Geometries::geometry6PPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 6;
	allReflections.resize(5);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1
	allReflections[0][3] = 4;
	allReflections[0][4] = 3;
	allReflections[0][2] = 5;
	allReflections[0][5] = 2;

	//P-2
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;
	allReflections[1][1] = 3;
	allReflections[1][3] = 1;

	//P-3
	allReflections[2][1] = 5;
	allReflections[2][5] = 1;
	allReflections[2][2] = 4;
	allReflections[2][4] = 2;

	//P-4
	allReflections[3][1] = 2;
	allReflections[3][2] = 1;
	allReflections[3][3] = 5;
	allReflections[3][5] = 3;

	//P-5
	allReflections[4][2] = 3;
	allReflections[4][3] = 2;
	allReflections[4][1] = 4;
	allReflections[4][4] = 1;

}

string Geometries::geometry6PPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C5-1-p";
			break;
		case 1:
			return "C5-1-pp";
			break;
		case 2:
			return "C5-1-mm";
			break;
		case 3:
			return "C5-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}



void Geometries::geometry7COCotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 7;
	allReflections.resize(3);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][1] = 6;
	allReflections[0][6] = 1;
	allReflections[0][0] = 3;
	allReflections[0][3] = 0;

	//P-2
	allReflections[1][1] = 4;
	allReflections[1][4] = 1;
	allReflections[1][2] = 3;
	allReflections[1][3] = 2;

	//P-3
	allReflections[2][4] = 6;
	allReflections[2][6] = 4;
	allReflections[2][0] = 2;
	allReflections[2][2] = 0;
}

string Geometries::geometry7COCSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
			break;

		case 1:
			return "C3-1-p";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}

void Geometries::geometry7CTPRotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 7;
	allReflections.resize(2);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;
	allReflections[0][3] = 4;
	allReflections[0][4] = 3;
	allReflections[0][5] = 6;
	allReflections[0][6] = 5;

	//P-2
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;
	allReflections[1][1] = 3;
	allReflections[1][3] = 1;
}

string Geometries::geometry7CTPRSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry7PBPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 7;
	allReflections.resize(10);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;

	//P-2
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;
	allReflections[1][6] = 5;
	allReflections[1][5] = 6;

	//P-3
	allReflections[2][2] = 6;
	allReflections[2][6] = 2;
	allReflections[2][3] = 5;
	allReflections[2][5] = 3;

	//P-4
	allReflections[3][2] = 3;
	allReflections[3][3] = 2;
	allReflections[3][6] = 4;
	allReflections[3][4] = 6;

	//P-5
	allReflections[4][3] = 4;
	allReflections[4][4] = 3;
	allReflections[4][2] = 5;
	allReflections[4][5] = 2;

	//P-6
	allReflections[5][4] = 5;
	allReflections[5][5] = 4;
	allReflections[5][3] = 6;
	allReflections[5][6] = 3;

	//S5-1-p
	allReflections[6][0] = 1;
	allReflections[6][1] = 0;
	allReflections[6][2] = 3;
	allReflections[6][3] = 4;
	allReflections[6][4] = 5;
	allReflections[6][5] = 6;
	allReflections[6][6] = 2;

	//S5-1-pp
	allReflections[7][0] = 1;
	allReflections[7][1] = 0;
	allReflections[7][2] = 4;
	allReflections[7][4] = 6;
	allReflections[7][6] = 3;
	allReflections[7][3] = 5;
	allReflections[7][5] = 2;

	//S5-1-m
	allReflections[8][0] = 1;
	allReflections[8][1] = 0;
	allReflections[8][2] = 6;
	allReflections[8][6] = 5;
	allReflections[8][5] = 4;
	allReflections[8][4] = 3;
	allReflections[8][3] = 2;

	//S5-1-mm
	allReflections[9][0] = 1;
	allReflections[9][1] = 0;
	allReflections[9][2] = 5;
	allReflections[9][5] = 3;
	allReflections[9][3] = 6;
	allReflections[9][6] = 4;
	allReflections[9][4] = 2;
}

string Geometries::geometry7PBPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C5-1-p";
		case 1:
			return "C5-1-pp";
		case 2:
			return "C5-1-m";
		case 3:
			return "C5-1-mm";
		case 4:
			return "C2-1 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )";
		case 5:
			return "C2-2 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )";
		case 6:
			return "C2-3 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )";
		case 7:
			return "C2-4 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )";
		case 8:
			return "C2-5 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C5-1-p C5-1-pp C5-1-m C5-1-mm )";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "S5-1-p";
		case 7:
			return "S5-1-pp";
		case 8:
			return "S5-1-m";
		case 9:
			return "S5-1-mm";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry7HPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 7;
	allReflections.resize(6);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][2] = 6;
	allReflections[0][6] = 2;
	allReflections[0][3] = 5;
	allReflections[0][5] = 3;

	//P-2
	allReflections[1][1] = 2;
	allReflections[1][2] = 1;
	allReflections[1][3] = 6;
	allReflections[1][6] = 3;
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;

	//P-3
	allReflections[2][1] = 3;
	allReflections[2][3] = 1;
	allReflections[2][4] = 6;
	allReflections[2][6] = 4;

	//P-4
	allReflections[3][2] = 3;
	allReflections[3][3] = 2;
	allReflections[3][1] = 4;
	allReflections[3][4] = 1;
	allReflections[3][6] = 5;
	allReflections[3][5] = 6;

	//P-5
	allReflections[4][2] = 4;
	allReflections[4][4] = 2;
	allReflections[4][1] = 5;
	allReflections[4][5] = 1;

	//P-6
	allReflections[5][4] = 3;
	allReflections[5][3] = 4;
	allReflections[5][2] = 5;
	allReflections[5][5] = 2;
	allReflections[5][1] = 6;
	allReflections[5][6] = 1;

}

string Geometries::geometry7HPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C6-1-p ";
		case 1:
			return "C3-1-p ";
		case 2:
			return "C2-1 ";
		case 3:
			return "C3-1-m ";
		case 4:
			return "C6-1-m ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry7HPotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 7;
	allReflections.resize(14);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1

	//P-2
	allReflections[1][4] = 3;
	allReflections[1][3] = 4;
	allReflections[1][2] = 5;
	allReflections[1][5] = 2;
	allReflections[1][1] = 6;
	allReflections[1][6] = 1;

	//P-3
	allReflections[2][4] = 5;
	allReflections[2][5] = 4;
	allReflections[2][3] = 6;
	allReflections[2][6] = 3;
	allReflections[2][0] = 2;
	allReflections[2][2] = 0;

	//P-4
	allReflections[3][6] = 5;
	allReflections[3][5] = 6;
	allReflections[3][0] = 4;
	allReflections[3][4] = 0;
	allReflections[3][1] = 3;
	allReflections[3][3] = 1;

	//P-5
	allReflections[4][0] = 6;
	allReflections[4][6] = 0;
	allReflections[4][1] = 5;
	allReflections[4][5] = 1;
	allReflections[4][2] = 4;
	allReflections[4][4] = 2;

	//P-6
	allReflections[5][0] = 1;
	allReflections[5][1] = 0;
	allReflections[5][2] = 6;
	allReflections[5][6] = 2;
	allReflections[5][3] = 5;
	allReflections[5][5] = 3;

	//P-7
	allReflections[6][2] = 1;
	allReflections[6][1] = 2;
	allReflections[6][0] = 3;
	allReflections[6][3] = 0;
	allReflections[6][4] = 6;
	allReflections[6][6] = 4;

	//P-8
	allReflections[7][2] = 3;
	allReflections[7][3] = 2;
	allReflections[7][1] = 4;
	allReflections[7][4] = 1;
	allReflections[7][0] = 5;
	allReflections[7][5] = 0;

	//S7-1-p
	allReflections[8][0] = 1;
	allReflections[8][1] = 2;
	allReflections[8][2] = 3;
	allReflections[8][3] = 4;
	allReflections[8][4] = 5;
	allReflections[8][5] = 6;
	allReflections[8][6] = 0;

	//S7-1-pp
	allReflections[9][0] = 2;
	allReflections[9][2] = 4;
	allReflections[9][4] = 6;
	allReflections[9][6] = 1;
	allReflections[9][1] = 3;
	allReflections[9][3] = 5;
	allReflections[9][5] = 0;

	//S7-1-ppp
	allReflections[10][0] = 3;
	allReflections[10][3] = 6;
	allReflections[10][6] = 2;
	allReflections[10][2] = 5;
	allReflections[10][5] = 1;
	allReflections[10][1] = 4;
	allReflections[10][4] = 0;

	//S7-1-m
	allReflections[11][0] = 6;
	allReflections[11][6] = 5;
	allReflections[11][5] = 4;
	allReflections[11][4] = 3;
	allReflections[11][3] = 2;
	allReflections[11][2] = 1;
	allReflections[11][1] = 0;

	//S7-1-mm
	allReflections[12][0] = 5;
	allReflections[12][5] = 3;
	allReflections[12][3] = 1;
	allReflections[12][1] = 6;
	allReflections[12][6] = 4;
	allReflections[12][4] = 2;
	allReflections[12][2] = 0;

	//S7-1-mmm
	allReflections[13][0] = 4;
	allReflections[13][4] = 1;
	allReflections[13][1] = 5;
	allReflections[13][5] = 2;
	allReflections[13][2] = 6;
	allReflections[13][6] = 3;
	allReflections[13][3] = 0;

}

string Geometries::geometry7HPSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{ 
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C7-1-p ";
		case 1:
			return "C7-1-pp ";
		case 2:
			return "C7-1-ppp ";
		case 3:
			return "C7-1-mmm ";
		case 4:
			return "C7-1-mm ";
		case 5:
			return "C7-1-m ";
		case 6:
			return "C2-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 7:
			return "C2-2 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 8:
			return "C2-3 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 9:
			return "C2-4 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 10:
			return "C2-5 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 11:
			return "C2-6 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 12:
			return "C2-7 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-mmm C7-1-mm C7-1-m ) ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "P-7 ";
		case 7:
			return "P-8 ";
		case 8:
			return "S7-1-p";
		case 9:
			return "S7-1-pp";
		case 10:
			return "S7-1-ppp";
		case 11:
			return "S7-1-m";
		case 12:
			return "S7-1-mm";
		case 13:
			return "S7-1-mmm";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry7JETPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 7;
	allReflections.resize(3);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][0] = 4;
	allReflections[0][4] = 0;
	allReflections[0][1] = 5;
	allReflections[0][5] = 1;

	//P-2
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;
	allReflections[1][3] = 5;
	allReflections[1][5] = 3;

	//P-3
	allReflections[2][0] = 2;
	allReflections[2][2] = 0;
	allReflections[2][1] = 3;
	allReflections[2][3] = 1;

}

string Geometries::geometry7JETPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p ";
		case 1:
			return "C3-1-m ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry8HBPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(12);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;

	//P-2
	allReflections[1][2] = 7;
	allReflections[1][7] = 2;
	allReflections[1][3] = 6;
	allReflections[1][6] = 3;
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;

	//P-3
	allReflections[2][4] = 6;
	allReflections[2][6] = 4;
	allReflections[2][3] = 7;
	allReflections[2][7] = 3;

	//P-4
	allReflections[3][6] = 5;
	allReflections[3][5] = 6;
	allReflections[3][4] = 7;
	allReflections[3][7] = 4;
	allReflections[3][2] = 3;
	allReflections[3][3] = 2;

	//P-5
	allReflections[4][5] = 7;
	allReflections[4][7] = 5;
	allReflections[4][2] = 4;
	allReflections[4][4] = 2;

	//P-6
	allReflections[5][3] = 4;
	allReflections[5][4] = 3;
	allReflections[5][2] = 5;
	allReflections[5][5] = 2;
	allReflections[5][6] = 7;
	allReflections[5][7] = 6;

	//P-7
	allReflections[6][2] = 6;
	allReflections[6][6] = 2;
	allReflections[6][3] = 5;
	allReflections[6][5] = 3;

	//Inv
	allReflections[7][0] = 1;
	allReflections[7][1] = 0;
	allReflections[7][4] = 7;
	allReflections[7][7] = 4;
	allReflections[7][3] = 6;
	allReflections[7][6] = 3;
	allReflections[7][2] = 5;
	allReflections[7][5] = 2;

	//S6-1-p
	allReflections[8][0] = 1;
	allReflections[8][1] = 0;
	allReflections[8][2] = 3;
	allReflections[8][3] = 4;
	allReflections[8][4] = 5;
	allReflections[8][5] = 6;
	allReflections[8][6] = 7;
	allReflections[8][7] = 2;

	//S6-1-m
	allReflections[9][0] = 1;
	allReflections[9][1] = 0;
	allReflections[9][2] = 7;
	allReflections[9][7] = 6;
	allReflections[9][6] = 5;
	allReflections[9][5] = 4;
	allReflections[9][4] = 3;
	allReflections[9][3] = 2;

	//S3-1-p
	allReflections[10][0] = 1;
	allReflections[10][1] = 0;
	allReflections[10][2] = 4;
	allReflections[10][4] = 6;
	allReflections[10][6] = 2;
	allReflections[10][3] = 5;
	allReflections[10][5] = 7;
	allReflections[10][7] = 3;

	//S3-1-m
	allReflections[11][0] = 1;
	allReflections[11][1] = 0;
	allReflections[11][2] = 6;
	allReflections[11][6] = 4;
	allReflections[11][4] = 2;
	allReflections[11][3] = 7;
	allReflections[11][7] = 5;
	allReflections[11][5] = 3;

}

string Geometries::geometry8HBPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C6-1-p";
		case 1:
			return "C3-1-p";
		case 2:
			return "C2-1 ";
		case 3:
			return "C3-1-m";
		case 4:
			return "C6-1-m";
		case 5:
			return "C2-2 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		case 6:
			return "C2-3 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		case 7:
			return "C2-4 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		case 8:
			return "C2-5 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		case 9:
			return "C2-6 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		case 10:
			return "C2-7 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C6-1-p C3-1-p C2-1 C6-1-m C3-1-m )";
		case 1:
			return "P-2 ( C2-3 )";
		case 2:
			return "P-3 ( C2-7 )";
		case 3:
			return "P-4 ( C2-4 )";
		case 4:
			return "P-5 ( C2-5 )";
		case 5:
			return "P-6 ( C2-2 )";
		case 6:
			return "P-7 ( C2-6 )";
		case 7:
			return "Inv";
		case 8:
			return "S6-1-p";
		case 9:
			return "S6-1-m";
		case 10:
			return "S3-1-p";
		case 11:
			return "S3-1-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}




void Geometries::geometry8BTPRotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(2);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;
	allReflections[0][2] = 3;
	allReflections[0][3] = 2;
	allReflections[0][4] = 5;
	allReflections[0][5] = 4;

	//P-2
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;
	allReflections[1][6] = 7;
	allReflections[1][7] = 6;
	allReflections[1][3] = 5;
	allReflections[1][5] = 3;
}


string Geometries::geometry8BTPRSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry8TDDotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(4);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][0] = 2;
	allReflections[0][2] = 0;
	allReflections[0][4] = 6;
	allReflections[0][6] = 4;

	//P-2
	allReflections[1][5] = 7;
	allReflections[1][7] = 5;
	allReflections[1][1] = 3;
	allReflections[1][3] = 1;

	//S4-1-p
	allReflections[2][0] = 7;
	allReflections[2][7] = 2;
	allReflections[2][2] = 5;
	allReflections[2][5] = 0;
	allReflections[2][3] = 6;
	allReflections[2][6] = 1;
	allReflections[2][1] = 4;
	allReflections[2][4] = 3;

	//S4-1-m
	allReflections[3][0] = 5;
	allReflections[3][5] = 2;
	allReflections[3][2] = 7;
	allReflections[3][7] = 0;
	allReflections[3][1] = 6;
	allReflections[3][6] = 3;
	allReflections[3][3] = 4;
	allReflections[3][4] = 1;

}


string Geometries::geometry8TDDSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ( C2-2 C2-3 )";
		case 1:
			return "C2-2 ( C2-1 C2-3 )";
		case 2:
			return "C2-3 ( C2-1 C2-2 )";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "S4-1-p";
		case 3:
			return "S4-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}

void Geometries::geometry8CUotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(24);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1 ( C4-1-p C2-1 C4-1-m )
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;
	allReflections[0][3] = 5;
	allReflections[0][5] = 3;
	allReflections[0][6] = 7;
	allReflections[0][7] = 6;
	allReflections[0][2] = 4;
	allReflections[0][4] = 2;

	//P-2 (C4-2-p C2-2 C4-2-m )
	allReflections[1][0] = 2;
	allReflections[1][2] = 0;
	allReflections[1][3] = 6;
	allReflections[1][6] = 3;
	allReflections[1][5] = 7;
	allReflections[1][7] = 5;
	allReflections[1][1] = 4;
	allReflections[1][4] = 1;

	//P-3 (C4-3-p C2-3 C4-3-m )
	allReflections[2][0] = 3;
	allReflections[2][3] = 0;
	allReflections[2][2] = 6;
	allReflections[2][6] = 2;
	allReflections[2][4] = 7;
	allReflections[2][7] = 4;
	allReflections[2][1] = 5;
	allReflections[2][5] = 1;

	//P4 ( C2-4 )
	allReflections[3][0] = 6;
	allReflections[3][6] = 0;
	allReflections[3][1] = 7;
	allReflections[3][7] = 1;

	//P-5 ( C2-5 )
	allReflections[4][0] = 5;
	allReflections[4][5] = 0;
	allReflections[4][2] = 7;
	allReflections[4][7] = 2;

	//P-6 ( C2-6 )
	allReflections[5][0] = 4;
	allReflections[5][4] = 0;
	allReflections[5][3] = 7;
	allReflections[5][7] = 3;

	//P-7 ( C2-7 )
	allReflections[6][2] = 3;
	allReflections[6][3] = 2;
	allReflections[6][4] = 5;
	allReflections[6][5] = 4;

	//P-8 ( C2-8 )
	allReflections[7][1] = 3;
	allReflections[7][3] = 1;
	allReflections[7][4] = 6;
	allReflections[7][6] = 4;

	//P-9 ( C2-9 )
	allReflections[8][1] = 2;
	allReflections[8][2] = 1;
	allReflections[8][5] = 6;
	allReflections[8][6] = 5;

	//Inv
	allReflections[9][0] = 7;
	allReflections[9][7] = 0;
	allReflections[9][3] = 4;
	allReflections[9][4] = 3;
	allReflections[9][1] = 6;
	allReflections[9][6] = 1;
	allReflections[9][2] = 5;
	allReflections[9][5] = 2;

	//S4-1-p
	allReflections[10][0] = 4;
	allReflections[10][4] = 6;
	allReflections[10][6] = 5;
	allReflections[10][5] = 0;
	allReflections[10][2] = 7;
	allReflections[10][7] = 3;
	allReflections[10][3] = 1;
	allReflections[10][1] = 2;

	//S4-1-m
	allReflections[11][0] = 5;
	allReflections[11][5] = 6;
	allReflections[11][6] = 4;
	allReflections[11][4] = 0;
	allReflections[11][2] = 1;
	allReflections[11][1] = 3;
	allReflections[11][3] = 7;
	allReflections[11][7] = 2;

	//S4-2-p
	allReflections[12][0] = 6;
	allReflections[12][6] = 5;
	allReflections[12][5] = 4;
	allReflections[12][4] = 0;
	allReflections[12][3] = 7;
	allReflections[12][7] = 1;
	allReflections[12][1] = 2;
	allReflections[12][2] = 3;

	//S4-2-m
	allReflections[13][0] = 4;
	allReflections[13][4] = 5;
	allReflections[13][5] = 6;
	allReflections[13][6] = 0;
	allReflections[13][3] = 2;
	allReflections[13][2] = 1;
	allReflections[13][1] = 7;
	allReflections[13][7] = 3;

	//S4-3-p
	allReflections[14][0] = 5;
	allReflections[14][5] = 4;
	allReflections[14][4] = 6;
	allReflections[14][6] = 0;
	allReflections[14][1] = 7;
	allReflections[14][7] = 2;
	allReflections[14][2] = 3;
	allReflections[14][3] = 1;

	//S4-3-m
	allReflections[15][0] = 6;
	allReflections[15][6] = 4;
	allReflections[15][4] = 5;
	allReflections[15][5] = 0;
	allReflections[15][1] = 3;
	allReflections[15][3] = 2;
	allReflections[15][2] = 7;
	allReflections[15][7] = 1;

	//S6-1-p
	allReflections[16][0] = 7;
	allReflections[16][7] = 0;
	allReflections[16][1] = 4;
	allReflections[16][4] = 2;
	allReflections[16][2] = 6;
	allReflections[16][6] = 3;
	allReflections[16][3] = 5;
	allReflections[16][5] = 1;

	//S6-1-m
	allReflections[17][0] = 7;
	allReflections[17][7] = 0;
	allReflections[17][1] = 5;
	allReflections[17][5] = 3;
	allReflections[17][3] = 6;
	allReflections[17][6] = 2;
	allReflections[17][2] = 4;
	allReflections[17][4] = 1;

	//S6-2-p
	allReflections[18][1] = 6;
	allReflections[18][6] = 1;
	allReflections[18][0] = 3;
	allReflections[18][3] = 5;
	allReflections[18][5] = 7;
	allReflections[18][7] = 4;
	allReflections[18][4] = 2;
	allReflections[18][2] = 0;

	//S6-2-m
	allReflections[19][1] = 6;
	allReflections[19][6] = 1;
	allReflections[19][0] = 2;
	allReflections[19][2] = 4;
	allReflections[19][4] = 7;
	allReflections[19][7] = 5;
	allReflections[19][5] = 3;
	allReflections[19][3] = 0;

	//S6-3-p
	allReflections[20][2] = 5;
	allReflections[20][5] = 2;
	allReflections[20][0] = 1;
	allReflections[20][1] = 4;
	allReflections[20][4] = 7;
	allReflections[20][7] = 6;
	allReflections[20][6] = 3;
	allReflections[20][3] = 0;

	//S6-3-m
	allReflections[21][2] = 5;
	allReflections[21][5] = 2;
	allReflections[21][0] = 3;
	allReflections[21][3] = 6;
	allReflections[21][6] = 7;
	allReflections[21][7] = 4;
	allReflections[21][4] = 1;
	allReflections[21][1] = 0;

	//S6-4-p
	allReflections[22][3] = 4;
	allReflections[22][4] = 3;
	allReflections[22][0] = 1;
	allReflections[22][1] = 5;
	allReflections[22][5] = 7;
	allReflections[22][7] = 6;
	allReflections[22][6] = 2;
	allReflections[22][2] = 0;

	//S6-4-m
	allReflections[23][3] = 4;
	allReflections[23][4] = 3;
	allReflections[23][0] = 2;
	allReflections[23][2] = 6;
	allReflections[23][6] = 7;
	allReflections[23][7] = 5;
	allReflections[23][5] = 1;
	allReflections[23][1] = 0;

}


string Geometries::geometry8CUSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C4-1-p";
		case 1:
			return "C2-1 ( C4-2-p C2-2 C4-2-m C4-3-p C2-3 C4-3-m C2-4 C2-7 )";
		case 2:
			return "C4-1-m";
		case 3:
			return "C4-2-p";
		case 4:
			return "C2-2 ( C4-1-p C2-1 C4-1-m C4-3-p C2-3 C4-3-m C2-5 C2-8 )";
		case 5:
			return "C4-2-m";
		case 6:
			return "C4-3-p";
		case 7:
			return "C2-3 ( C4-1-p C2-1 C4-1-m C4-2-p C2-2 C4-2-m C2-6 C2-9 )";
		case 8:
			return "C4-3-m";
		case 9:
			return "C3-1-p";
		case 10:
			return "C3-1-m";
		case 11:
			return "C3-2-p";
		case 12:
			return "C3-2-m";
		case 13:
			return "C3-3-p";
		case 14:
			return "C3-3-m";
		case 15:
			return "C3-4-p";
		case 16:
			return "C3-4-m";
		case 17:
			return "C2-4 ( C3-3-p C3-3-m C3-4-p C3-4-m C4-1-p C2-1 C4-1-m )";
		case 18:
			return "C2-5 ( C3-2-p C3-2-m C3-4-p C4-4-m C4-2-p C2-2 C4-2-m )";
		case 19:
			return "C2-6 ( C3-2-p C3-2-m C3-3-p C3-3-m C4-3-p C2-3 C4-3-m )";
		case 20:
			return "C2-7 ( C3-1-p C3-1-m C3-2-p C3-2-m C4-1-p C2-1 C4-1-m )";
		case 21:
			return "C2-8 ( C3-1-p C3-1-m C3-3-p C3-3-m C4-2-p C2-2 C4-2-m )";
		case 22:
			return "C2-9 ( C3-1-p C3-1-m C3-4-p C3-4-m C4-3-p C2-3 C4-3-m )";

		default:
			cout << "rotation on Geometries::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C4-1-p C2-1 C4-1-m )";
		case 1:
			return "P-2 ( C4-2-p C2-2 C4-2-m )";
		case 2:
			return "P-3 ( C4-3-p C2-3 C4-3-m )";
		case 3:
			return "P-4 ( C2-4 )";
		case 4:
			return "P-5 ( C2-5 )";
		case 5:
			return "P-6 ( C2-6 )";
		case 6:
			return "P-7 ( C2-7 )";
		case 7:
			return "P-8 ( C2-8 )";
		case 8:
			return "P-9 ( C2-9 )";
		case 9:
			return "Inv";
		case 10:
			return "S4-1-p";
		case 11:
			return "S4-1-m";
		case 12:
			return "S4-2-p";
		case 13:
			return "S4-2-m";
		case 14:
			return "S4-3-p";
		case 15:
			return "S4-3-m";
		case 16:
			return "S6-1-p";
		case 17:
			return "S6-1-m";
		case 18:
			return "S6-2-p";
		case 19:
			return "S6-2-m";
		case 20:
			return "S6-3-p";
		case 21:
			return "S6-3-m";
		case 22:
			return "S6-4-p";
		case 23:
			return "S6-4-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry8SAPRotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(8);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1
	allReflections[0][3] = 4;
	allReflections[0][4] = 3;
	allReflections[0][6] = 7;
	allReflections[0][7] = 6;
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;

	//P-2
	allReflections[1][1] = 5;
	allReflections[1][5] = 1;
	allReflections[1][4] = 6;
	allReflections[1][6] = 4;
	allReflections[1][0] = 2;
	allReflections[1][2] = 0;

	//P-3
	allReflections[2][0] = 5;
	allReflections[2][5] = 0;
	allReflections[2][3] = 6;
	allReflections[2][6] = 3;
	allReflections[2][4] = 7;
	allReflections[2][7] = 4;

	//P-4
	allReflections[3][0] = 1;
	allReflections[3][1] = 0;
	allReflections[3][2] = 5;
	allReflections[3][5] = 2;
	allReflections[3][3] = 7;
	allReflections[3][7] = 3;

	//S8-1-p
	allReflections[4][0] = 3;
	allReflections[4][3] = 2;
	allReflections[4][2] = 6;
	allReflections[4][6] = 5;
	allReflections[4][5] = 7;
	allReflections[4][7] = 1;
	allReflections[4][1] = 4;
	allReflections[4][4] = 0;

	//S8-1-pp
	allReflections[5][0] = 6;
	allReflections[5][6] = 1;
	allReflections[5][1] = 3;
	allReflections[5][3] = 5;
	allReflections[5][5] = 4;
	allReflections[5][4] = 2;
	allReflections[5][2] = 7;
	allReflections[5][7] = 0;

	//S8-1-m
	allReflections[6][0] = 4;
	allReflections[6][4] = 1;
	allReflections[6][1] = 7;
	allReflections[6][7] = 5;
	allReflections[6][5] = 6;
	allReflections[6][6] = 2;
	allReflections[6][2] = 3;
	allReflections[6][3] = 0;

	//S8-1-mm
	allReflections[7][0] = 7;
	allReflections[7][7] = 2;
	allReflections[7][2] = 4;
	allReflections[7][4] = 5;
	allReflections[7][5] = 3;
	allReflections[7][3] = 1;
	allReflections[7][1] = 6;
	allReflections[7][6] = 0;

}


string Geometries::geometry8SAPRSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C4-1-p";
		case 1:
			return "C2-1 ";
		case 2:
			return "C4-1-m";
		case 3:
			return "C2-2 ( C4-1-p C2-1 C4-1-m )";
		case 4:
			return "C2-3 ( C4-1-p C2-1 C4-1-m )";
		case 5:
			return "C2-4 ( C4-1-p C2-1 C4-1-m )";
		case 6:
			return "C2-5 ( C4-1-p C2-1 C4-1-m )";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "S8-1-p";
		case 5:
			return "S8-1-pp";
		case 6:
			return "S8-1-m";
		case 7:
			return "S8-1-mm";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry8ETBPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(6);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1 ( C3-1-p C3-1-m )
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;
	allReflections[0][2] = 3;
	allReflections[0][3] = 2;
	allReflections[0][4] = 5;
	allReflections[0][5] = 4;
	allReflections[0][6] = 7;
	allReflections[0][7] = 6;

	//P-2
	allReflections[1][3] = 5;
	allReflections[1][5] = 3;
	allReflections[1][2] = 4;
	allReflections[1][4] = 2;

	//P-3
	allReflections[2][1] = 5;
	allReflections[2][5] = 1;
	allReflections[2][0] = 4;
	allReflections[2][4] = 0;

	//P-4
	allReflections[3][1] = 3;
	allReflections[3][3] = 1;
	allReflections[3][0] = 2;
	allReflections[3][2] = 0;

	//S3-1-p
	allReflections[4][6] = 7;
	allReflections[4][7] = 6;
	allReflections[4][1] = 2;
	allReflections[4][2] = 5;
	allReflections[4][5] = 0;
	allReflections[4][0] = 3;
	allReflections[4][3] = 4;
	allReflections[4][4] = 1;

	//S3-1-m
	allReflections[5][7] = 6;
	allReflections[5][6] = 7;
	allReflections[5][1] = 4;
	allReflections[5][4] = 3;
	allReflections[5][3] = 0;
	allReflections[5][0] = 5;
	allReflections[5][5] = 2;
	allReflections[5][2] = 1;

}


string Geometries::geometry8ETBPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
		case 1:
			return "C3-1-m";
		case 2:
			return "C2-1 ( C3-1-p C3-1-m ) ";
		case 3:
			return "C2-2 ( C3-1-p C3-1-m ) ";
		case 4:
			return "C2-3 ( C3-1-p C3-1-m ) ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C3-1-p C3-1-m ) ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "S3-1-p ";
		case 5:
			return "S3-1-m ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry8HPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(7);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;
	allReflections[0][3] = 7;
	allReflections[0][7] = 3;
	allReflections[0][4] = 6;
	allReflections[0][6] = 4;

	//P-2
	allReflections[1][2] = 3;
	allReflections[1][3] = 2;
	allReflections[1][1] = 4;
	allReflections[1][4] = 1;
	allReflections[1][5] = 7;
	allReflections[1][7] = 5;

	//P-3
	allReflections[2][3] = 4;
	allReflections[2][4] = 3;
	allReflections[2][2] = 5;
	allReflections[2][5] = 2;
	allReflections[2][1] = 6;
	allReflections[2][6] = 1;

	//P-4
	allReflections[3][4] = 5;
	allReflections[3][5] = 4;
	allReflections[3][3] = 6;
	allReflections[3][6] = 3;
	allReflections[3][2] = 7;
	allReflections[3][7] = 2;

	//P-5
	allReflections[4][6] = 5;
	allReflections[4][5] = 6;
	allReflections[4][4] = 7;
	allReflections[4][7] = 4;
	allReflections[4][1] = 3;
	allReflections[4][3] = 1;

	//P-6
	allReflections[5][6] = 7;
	allReflections[5][7] = 6;
	allReflections[5][1] = 5;
	allReflections[5][5] = 1;
	allReflections[5][2] = 4;
	allReflections[5][4] = 2;

	//P-7
	allReflections[6][1] = 7;
	allReflections[6][7] = 1;
	allReflections[6][2] = 6;
	allReflections[6][6] = 2;
	allReflections[6][3] = 5;
	allReflections[6][5] = 3;

}


string Geometries::geometry8HPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C7-1-p";
		case 1:
			return "C7-1-pp";
		case 2:
			return "C7-1-ppp";
		case 3:
			return "C7-1-mmm";
		case 4:
			return "C7-1-mm";
		case 5:
			return "C7-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "P-7 ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry8OPotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(16);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1  ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m ) 

	//P-2 ( C2-6 )
	allReflections[1][1] = 7;
	allReflections[1][7] = 1;
	allReflections[1][2] = 6;
	allReflections[1][6] = 2;
	allReflections[1][3] = 5;
	allReflections[1][5] = 3;

	//P-3 ( C2-7 )
	allReflections[2][0] = 1;
	allReflections[2][1] = 0;
	allReflections[2][2] = 7;
	allReflections[2][7] = 2;
	allReflections[2][3] = 6;
	allReflections[2][6] = 3;
	allReflections[2][4] = 5;
	allReflections[2][5] = 4;

	//P-4 ( C2-8 )
	allReflections[3][0] = 2;
	allReflections[3][2] = 0;
	allReflections[3][3] = 7;
	allReflections[3][7] = 3;
	allReflections[3][4] = 6;
	allReflections[3][6] = 4;

	//P-5 ( C2-9 )
	allReflections[4][1] = 2;
	allReflections[4][2] = 1;
	allReflections[4][0] = 3;
	allReflections[4][3] = 0;
	allReflections[4][4] = 7;
	allReflections[4][7] = 4;
	allReflections[4][6] = 5;
	allReflections[4][5] = 6;

	//P-6 ( C2-2 )
	allReflections[5][1] = 3;
	allReflections[5][3] = 1;
	allReflections[5][0] = 4;
	allReflections[5][4] = 0;
	allReflections[5][5] = 7;
	allReflections[5][7] = 5;

	//P-7 ( C2-3 )
	allReflections[6][2] = 3;
	allReflections[6][3] = 2;
	allReflections[6][1] = 4;
	allReflections[6][4] = 1;
	allReflections[6][0] = 5;
	allReflections[6][5] = 0;
	allReflections[6][6] = 7;
	allReflections[6][7] = 6;

	//P-8 ( C2-4 )
	allReflections[7][2] = 4;
	allReflections[7][4] = 2;
	allReflections[7][1] = 5;
	allReflections[7][5] = 1;
	allReflections[7][0] = 6;
	allReflections[7][6] = 0;

	//P-9 ( C2-5 )
	allReflections[8][3] = 4;
	allReflections[8][4] = 3;
	allReflections[8][2] = 5;
	allReflections[8][5] = 2;
	allReflections[8][1] = 6;
	allReflections[8][6] = 1;
	allReflections[8][0] = 7;
	allReflections[8][7] = 0;

	//S8-1-p
	allReflections[9][0] = 1;
	allReflections[9][1] = 2;
	allReflections[9][2] = 3;
	allReflections[9][3] = 4;
	allReflections[9][4] = 5;
	allReflections[9][5] = 6;
	allReflections[9][6] = 7;
	allReflections[9][7] = 0;

	//S4-1-p
	allReflections[10][0] = 2;
	allReflections[10][2] = 4;
	allReflections[10][4] = 6;
	allReflections[10][6] = 0;
	allReflections[10][1] = 3;
	allReflections[10][3] = 5;
	allReflections[10][5] = 7;
	allReflections[10][7] = 1;

	//S8-1-ppp
	allReflections[11][0] = 3;
	allReflections[11][3] = 6;
	allReflections[11][6] = 1;
	allReflections[11][1] = 4;
	allReflections[11][4] = 7;
	allReflections[11][7] = 2;
	allReflections[11][2] = 5;
	allReflections[11][5] = 0;

	//S8-1-mmm
	allReflections[12][0] = 5;
	allReflections[12][5] = 2;
	allReflections[12][2] = 7;
	allReflections[12][7] = 4;
	allReflections[12][4] = 1;
	allReflections[12][1] = 6;
	allReflections[12][6] = 3;
	allReflections[12][3] = 0;

	//S4-1-m
	allReflections[13][0] = 6;
	allReflections[13][6] = 4;
	allReflections[13][4] = 2;
	allReflections[13][2] = 0;
	allReflections[13][1] = 7;
	allReflections[13][7] = 5;
	allReflections[13][5] = 3;
	allReflections[13][3] = 1;

	//S8-1-m
	allReflections[14][0] = 7;
	allReflections[14][7] = 6;
	allReflections[14][6] = 5;
	allReflections[14][5] = 4;
	allReflections[14][4] = 3;
	allReflections[14][3] = 2;
	allReflections[14][2] = 1;
	allReflections[14][1] = 0;

	//Inv
	allReflections[15][0] = 4;
	allReflections[15][4] = 0;
	allReflections[15][1] = 5;
	allReflections[15][5] = 1;
	allReflections[15][2] = 6;
	allReflections[15][6] = 2;
	allReflections[15][3] = 7;
	allReflections[15][7] = 3;


}


string Geometries::geometry8OPSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C8-1-p";
		case 1:
			return "C4-1-p";
		case 2:
			return "C8-1-ppp";
		case 3:
			return "C2-1 ( C2-2 C2-3 C2-4 C2-5 C2-6 C2-7 C2-8 C2-9 ) ";
		case 4:
			return "C8-1-mmm";
		case 5:
			return "C4-1-m";
		case 6:
			return "C8-1-m";
		case 7:
			return "C2-2 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-6 ) ";
		case 8:
			return "C2-3 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-7 ) ";
		case 9:
			return "C2-4 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-8 ) ";
		case 10:
			return "C2-5 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-9 ) ";
		case 11:
			return "C2-6 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-2 ) ";
		case 12:
			return "C2-7 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-3 ) ";
		case 13:
			return "C2-8 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-4 ) ";
		case 14:
			return "C2-9 ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m C2-5 ) ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1  ( C8-1-p C4-1-p C8-1-ppp C2-1 C8-1-mmm C4-1-m C8-1-m ) ";
		case 1:
			return "P-2 ( C2-6 ) ";
		case 2:
			return "P-3 ( C2-7 ) ";
		case 3:
			return "P-4 ( C2-8 ) ";
		case 4:
			return "P-5 ( C2-9 ) ";
		case 5:
			return "P-6 ( C2-2 ) ";
		case 6:
			return "P-7 ( C2-3 ) ";
		case 7:
			return "P-8 ( C2-4 ) ";
		case 8:
			return "P-9 ( C2-5 ) ";
		case 9:
			return "S8-1-p ";
		case 10:
			return "S4-1-p ";
		case 11:
			return "S8-1-ppp ";
		case 12:
			return "S8-1-mmm ";
		case 13:
			return "S4-1-m ";
		case 14:
			return "S8-1-m ";
		case 15:
			return "Inv ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry8JGBFotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 8;
	allReflections.resize(4);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;
	allReflections[0][4] = 3;
	allReflections[0][3] = 4;
	allReflections[0][2] = 5;
	allReflections[0][5] = 2;

	//P-2
	allReflections[1][2] = 3;
	allReflections[1][3] = 2;
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;
	allReflections[1][6] = 7;
	allReflections[1][7] = 6;

	//S4-1-p
	allReflections[2][1] = 6;
	allReflections[2][6] = 0;
	allReflections[2][0] = 7;
	allReflections[2][7] = 1;
	allReflections[2][3] = 4;
	allReflections[2][4] = 5;
	allReflections[2][5] = 2;
	allReflections[2][2] = 3;

	//S4-1-m
	allReflections[3][1] = 7;
	allReflections[3][7] = 0;
	allReflections[3][0] = 6;
	allReflections[3][6] = 1;
	allReflections[3][4] = 3;
	allReflections[3][3] = 2;
	allReflections[3][2] = 5;
	allReflections[3][5] = 4;

}


string Geometries::geometry8JGBFSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ( C2-2 C2-3 ) ";
		case 1:
			return "C2-2 ( C2-1 C2-3 ) ";
		case 2:
			return "C2-3 ( C2-1 C2-2 ) ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "S4-1-p";
		case 3:
			return "S4-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}





void Geometries::geometry9TCTPRotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(6);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1 ( C3-1-p C3-1-m )
	allReflections[0][2] = 3;
	allReflections[0][3] = 2;
	allReflections[0][7] = 8;
	allReflections[0][8] = 7;
	allReflections[0][1] = 4;
	allReflections[0][4] = 1;

	//P-2
	allReflections[1][3] = 4;
	allReflections[1][4] = 3;
	allReflections[1][1] = 2;
	allReflections[1][2] = 1;
	allReflections[1][6] = 5;
	allReflections[1][5] = 6;

	//P-3
	allReflections[2][3] = 8;
	allReflections[2][8] = 3;
	allReflections[2][2] = 7;
	allReflections[2][7] = 2;
	allReflections[2][0] = 5;
	allReflections[2][5] = 0;

	//P-4
	allReflections[3][4] = 8;
	allReflections[3][8] = 4;
	allReflections[3][1] = 7;
	allReflections[3][7] = 1;
	allReflections[3][0] = 6;
	allReflections[3][6] = 0;

	//S3-1-p
	allReflections[4][1] = 3;
	allReflections[4][3] = 7;
	allReflections[4][7] = 4;
	allReflections[4][4] = 2;
	allReflections[4][2] = 8;
	allReflections[4][8] = 1;
	allReflections[4][0] = 6;
	allReflections[4][6] = 5;
	allReflections[4][5] = 0;

	//S3-1-m
	allReflections[5][1] = 8;
	allReflections[5][8] = 2;
	allReflections[5][2] = 4;
	allReflections[5][4] = 7;
	allReflections[5][7] = 3;
	allReflections[5][3] = 1;
	allReflections[5][0] = 5;
	allReflections[5][5] = 6;
	allReflections[5][6] = 0;

}


string Geometries::geometry9TCTPRSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ( C3-1-p C3-1-m )";
		case 1:
			return "C2-2 ( C3-1-p C3-1-m )";
		case 2:
			return "C2-3 ( C3-1-p C3-1-m )";
		case 3:
			return "C3-1-p";
		case 4:
			return "C3-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C3-1-p C3-1-m )";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "S3-1-p";
		case 5:
			return "S3-1-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry9CSAPRotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(4);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][5] = 6;
	allReflections[0][6] = 5;
	allReflections[0][1] = 3;
	allReflections[0][3] = 1;
	allReflections[0][8] = 7;
	allReflections[0][7] = 8;

	//P-2
	allReflections[1][1] = 2;
	allReflections[1][2] = 1;
	allReflections[1][8] = 6;
	allReflections[1][6] = 8;
	allReflections[1][4] = 3;
	allReflections[1][3] = 4;

	//P-3
	allReflections[2][8] = 5;
	allReflections[2][5] = 8;
	allReflections[2][4] = 2;
	allReflections[2][2] = 4;
	allReflections[2][7] = 6;
	allReflections[2][6] = 7;

	//P-4
	allReflections[3][1] = 4;
	allReflections[3][4] = 1;
	allReflections[3][7] = 5;
	allReflections[3][5] = 7;
	allReflections[3][3] = 2;
	allReflections[3][2] = 3;


}


string Geometries::geometry9CSAPRSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C4-1-p";
			break;
		case 1:
			return "C2-1 ";
			break;
		case 2:
			return "C4-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}




void Geometries::geometry9MFFotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(1);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][1] = 4;
	allReflections[0][4] = 1;
	allReflections[0][2] = 3;
	allReflections[0][3] = 2;
	allReflections[0][5] = 6;
	allReflections[0][6] = 5;

}


string Geometries::geometry9MFFSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		cout << "rotation on CauchyIndex::rotationString not found" << endl;
		exit(1);
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry9CCUotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(4);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][2] = 3;
	allReflections[0][3] = 2;
	allReflections[0][6] = 7;
	allReflections[0][7] = 6;
	allReflections[0][0] = 1;
	allReflections[0][1] = 0;
	allReflections[0][4] = 5;
	allReflections[0][5] = 4;

	//P-2
	allReflections[1][1] = 2;
	allReflections[1][2] = 1;
	allReflections[1][5] = 6;
	allReflections[1][6] = 5;

	//P-3
	allReflections[2][1] = 3;
	allReflections[2][3] = 1;
	allReflections[2][5] = 7;
	allReflections[2][7] = 5;
	allReflections[2][0] = 2;
	allReflections[2][2] = 0;
	allReflections[2][4] = 6;
	allReflections[2][6] = 4;

	//P-4
	allReflections[3][0] = 3;
	allReflections[3][3] = 0;
	allReflections[3][4] = 7;
	allReflections[3][7] = 4;

}


string Geometries::geometry9CCUSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C4-1-p";
			break;
		case 1:
			return "C2-1 ";
			break;
		case 2:
			return "C4-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry9HHotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(2);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][2] = 4;
	allReflections[0][4] = 2;
	allReflections[0][1] = 5;
	allReflections[0][5] = 1;

	//P-2
	allReflections[1][1] = 2;
	allReflections[1][2] = 1;
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;
	allReflections[1][0] = 3;
	allReflections[1][3] = 0;
	allReflections[1][8] = 7;
	allReflections[1][7] = 8;

}


string Geometries::geometry9HHSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C2-1 ";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry9OPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(8);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][1] = 7;
	allReflections[0][7] = 1;
	allReflections[0][2] = 6;
	allReflections[0][6] = 2;
	allReflections[0][3] = 5;
	allReflections[0][5] = 3;

	//P-2
	allReflections[1][1] = 8;
	allReflections[1][8] = 1;
	allReflections[1][2] = 7;
	allReflections[1][7] = 2;
	allReflections[1][3] = 6;
	allReflections[1][6] = 3;
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;

	//P-3
	allReflections[2][2] = 8;
	allReflections[2][8] = 2;
	allReflections[2][3] = 7;
	allReflections[2][7] = 3;
	allReflections[2][4] = 6;
	allReflections[2][6] = 4;

	//P-4
	allReflections[3][1] = 2;
	allReflections[3][2] = 1;
	allReflections[3][3] = 8;
	allReflections[3][8] = 3;
	allReflections[3][4] = 7;
	allReflections[3][7] = 4;
	allReflections[3][5] = 6;
	allReflections[3][6] = 5;

	//P-5
	allReflections[4][1] = 3;
	allReflections[4][3] = 1;
	allReflections[4][4] = 8;
	allReflections[4][8] = 4;
	allReflections[4][5] = 7;
	allReflections[4][7] = 5;

	//P-6
	allReflections[5][2] = 3;
	allReflections[5][3] = 2;
	allReflections[5][4] = 1;
	allReflections[5][1] = 4;
	allReflections[5][5] = 8;
	allReflections[5][8] = 5;
	allReflections[5][6] = 7;
	allReflections[5][7] = 6;

	//P-7
	allReflections[6][2] = 4;
	allReflections[6][4] = 2;
	allReflections[6][1] = 5;
	allReflections[6][5] = 1;
	allReflections[6][6] = 8;
	allReflections[6][8] = 6;

	//P-8
	allReflections[7][3] = 4;
	allReflections[7][4] = 3;
	allReflections[7][2] = 5;
	allReflections[7][5] = 2;
	allReflections[7][1] = 6;
	allReflections[7][6] = 1;
	allReflections[7][7] = 8;
	allReflections[7][8] = 7;

}


string Geometries::geometry9OPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C8-1-p";
			break;
		case 1:
			return "C4-1-p";
			break;
		case 2:
			return "C8-1-ppp";
			break;
		case 3:
			return "C2-1 ";
			break;
		case 4:
			return "C8-1-m";
			break;
		case 5:
			return "C4-1-m";
			break;
		case 6:
			return "C8-1-mmm";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "P-7 ";
		case 7:
			return "P-8 ";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry9EPotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(18);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )

	//P-2
	allReflections[1][3] = 4;
	allReflections[1][4] = 3;
	allReflections[1][2] = 5;
	allReflections[1][5] = 2;
	allReflections[1][1] = 6;
	allReflections[1][6] = 1;
	allReflections[1][0] = 7;
	allReflections[1][7] = 0;

	//P-3
	allReflections[2][2] = 3;
	allReflections[2][3] = 2;
	allReflections[2][1] = 4;
	allReflections[2][4] = 1;
	allReflections[2][0] = 5;
	allReflections[2][5] = 0;
	allReflections[2][8] = 6;
	allReflections[2][6] = 8;

	//P-4
	allReflections[3][1] = 2;
	allReflections[3][2] = 1;
	allReflections[3][0] = 3;
	allReflections[3][3] = 0;
	allReflections[3][8] = 4;
	allReflections[3][4] = 8;
	allReflections[3][7] = 5;
	allReflections[3][5] = 7;

	//P-5
	allReflections[4][0] = 1;
	allReflections[4][1] = 0;
	allReflections[4][8] = 2;
	allReflections[4][2] = 8;
	allReflections[4][7] = 3;
	allReflections[4][3] = 7;
	allReflections[4][6] = 4;
	allReflections[4][4] = 6;

	//P-6
	allReflections[5][0] = 8;
	allReflections[5][8] = 0;
	allReflections[5][1] = 7;
	allReflections[5][7] = 1;
	allReflections[5][2] = 6;
	allReflections[5][6] = 2;
	allReflections[5][3] = 5;
	allReflections[5][5] = 3;

	//P-7
	allReflections[6][8] = 7;
	allReflections[6][7] = 8;
	allReflections[6][0] = 6;
	allReflections[6][6] = 0;
	allReflections[6][1] = 5;
	allReflections[6][5] = 1;
	allReflections[6][2] = 4;
	allReflections[6][4] = 2;

	//P-8
	allReflections[7][6] = 7;
	allReflections[7][7] = 6;
	allReflections[7][5] = 8;
	allReflections[7][8] = 5;
	allReflections[7][4] = 0;
	allReflections[7][0] = 4;
	allReflections[7][1] = 3;
	allReflections[7][3] = 1;

	//P-9
	allReflections[8][6] = 5;
	allReflections[8][5] = 6;
	allReflections[8][7] = 4;
	allReflections[8][4] = 7;
	allReflections[8][8] = 3;
	allReflections[8][3] = 8;
	allReflections[8][0] = 2;
	allReflections[8][2] = 0;

	//P-10
	allReflections[9][4] = 5;
	allReflections[9][5] = 4;
	allReflections[9][3] = 6;
	allReflections[9][6] = 3;
	allReflections[9][2] = 7;
	allReflections[9][7] = 2;
	allReflections[9][1] = 8;
	allReflections[9][8] = 1;

	
	//S9-1-p
	allReflections[10][0] = 1;
	allReflections[10][1] = 2;
	allReflections[10][2] = 3;
	allReflections[10][3] = 4;
	allReflections[10][4] = 5;
	allReflections[10][5] = 6;
	allReflections[10][6] = 7;
	allReflections[10][7] = 8;
	allReflections[10][8] = 0;

	//S9-1-pp
	allReflections[11][0] = 2;
	allReflections[11][2] = 4;
	allReflections[11][4] = 6;
	allReflections[11][6] = 8;
	allReflections[11][8] = 1;
	allReflections[11][1] = 3;
	allReflections[11][3] = 5;
	allReflections[11][5] = 7;
	allReflections[11][7] = 0;

	//S3-1-p
	allReflections[12][0] = 3;
	allReflections[12][3] = 6;
	allReflections[12][6] = 0;
	allReflections[12][1] = 4;
	allReflections[12][4] = 7;
	allReflections[12][7] = 1;
	allReflections[12][2] = 5;
	allReflections[12][5] = 8;
	allReflections[12][8] = 2;

	//S9-1-ppp
	allReflections[13][0] = 4;
	allReflections[13][4] = 8;
	allReflections[13][8] = 3;
	allReflections[13][3] = 7;
	allReflections[13][7] = 2;
	allReflections[13][2] = 6;
	allReflections[13][6] = 1;
	allReflections[13][1] = 5;
	allReflections[13][5] = 0;

	//S9-1-mmm
	allReflections[14][0] = 5;
	allReflections[14][5] = 1;
	allReflections[14][1] = 6;
	allReflections[14][6] = 2;
	allReflections[14][2] = 7;
	allReflections[14][7] = 3;
	allReflections[14][3] = 8;
	allReflections[14][8] = 4;
	allReflections[14][4] = 0;

	//S3-1-m
	allReflections[15][0] = 6;
	allReflections[15][6] = 3;
	allReflections[15][3] = 0;
	allReflections[15][1] = 7;
	allReflections[15][7] = 4;
	allReflections[15][4] = 1;
	allReflections[15][2] = 8;
	allReflections[15][8] = 5;
	allReflections[15][5] = 2;

	//S9-1-mm
	allReflections[16][0] = 7;
	allReflections[16][7] = 5;
	allReflections[16][5] = 3;
	allReflections[16][3] = 1;
	allReflections[16][1] = 8;
	allReflections[16][8] = 6;
	allReflections[16][6] = 4;
	allReflections[16][4] = 2;
	allReflections[16][2] = 0;

	//S9-1-m
	allReflections[17][0] = 8;
	allReflections[17][8] = 7;
	allReflections[17][7] = 6;
	allReflections[17][6] = 5;
	allReflections[17][5] = 4;
	allReflections[17][4] = 3;
	allReflections[17][3] = 2;
	allReflections[17][2] = 1;
	allReflections[17][1] = 0;



}


string Geometries::geometry9EPSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C9-1-p";
			break;
		case 1:
			return "C9-1-pp";
			break;
		case 2:
			return "C3-1-p";
			break;
		case 3:
			return "C9-1-ppp";
			break;
		case 4:
			return "C9-1-mmm";
			break;
		case 5:
			return "C3-1-m";
			break;
		case 6:
			return "C9-1-mm";
			break;
		case 7:
			return "C9-1-m";
			break;
		case 8:
			return "C2-1 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 9:
			return "C2-2 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 10:
			return "C2-3 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 11:
			return "C2-4 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 12:
			return "C2-5 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 13:
			return "C2-6 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 14:
			return "C2-7 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 15:
			return "C2-8 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
			break;
		case 16:
			return "C2-9 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m )";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C9-1-p C9-1-pp C9-1-ppp C9-1-m C9-1-mm C9-1-mmm C3-1-p C3-1-m ) ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "P-7 ";
		case 7:
			return "P-8 ";
		case 8:
			return "P-9 ";
		case 9:
			return "P-10 ";
		case 10:
			return "S9-1-p";
		case 11:
			return "S9-1-pp";
		case 12:
			return "S3-1-p";
		case 13:
			return "S9-1-ppp";
		case 14:
			return "S9-1-mmm";
		case 15:
			return "S3-1-m";
		case 16:
			return "S9-1-mm";
		case 17:
			return "S9-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry9HBPYotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(14);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm )
	allReflections[0][0] = 8;
	allReflections[0][8] = 0;

	//P-2
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;
	allReflections[1][3] = 6;
	allReflections[1][6] = 3;
	allReflections[1][2] = 7;
	allReflections[1][7] = 2;

	//P-3
	allReflections[2][6] = 5;
	allReflections[2][5] = 6;
	allReflections[2][4] = 7;
	allReflections[2][7] = 4;
	allReflections[2][1] = 3;
	allReflections[2][3] = 1;

	//P-4
	allReflections[3][6] = 7;
	allReflections[3][7] = 6;
	allReflections[3][1] = 5;
	allReflections[3][5] = 1;
	allReflections[3][2] = 4;
	allReflections[3][4] = 2;

	//P-5
	allReflections[4][1] = 7;
	allReflections[4][7] = 1;
	allReflections[4][2] = 6;
	allReflections[4][6] = 2;
	allReflections[4][3] = 5;
	allReflections[4][5] = 3;

	//P-6
	allReflections[5][1] = 2;
	allReflections[5][2] = 1;
	allReflections[5][3] = 7;
	allReflections[5][7] = 3;
	allReflections[5][4] = 6;
	allReflections[5][6] = 4;

	//P-7
	allReflections[6][2] = 3;
	allReflections[6][3] = 2;
	allReflections[6][1] = 4;
	allReflections[6][4] = 1;
	allReflections[6][5] = 7;
	allReflections[6][7] = 5;

	//P-8
	allReflections[7][3] = 4;
	allReflections[7][4] = 3;
	allReflections[7][2] = 5;
	allReflections[7][5] = 2;
	allReflections[7][1] = 6;
	allReflections[7][6] = 1;


	//S7-1-p
	allReflections[8][0] = 8;
	allReflections[8][8] = 0;
	allReflections[8][1] = 7;
	allReflections[8][7] = 6;
	allReflections[8][6] = 5;
	allReflections[8][5] = 4;
	allReflections[8][4] = 3;
	allReflections[8][3] = 2;
	allReflections[8][2] = 1;

	//S7-1-pp
	allReflections[9][0] = 8;
	allReflections[9][8] = 0;
	allReflections[9][1] = 6;
	allReflections[9][6] = 4;
	allReflections[9][4] = 2;
	allReflections[9][2] = 7;
	allReflections[9][7] = 5;
	allReflections[9][5] = 3;
	allReflections[9][3] = 1;

	//S7-1-ppp
	allReflections[10][0] = 8;
	allReflections[10][8] = 0;
	allReflections[10][1] = 5;
	allReflections[10][5] = 2;
	allReflections[10][2] = 6;
	allReflections[10][6] = 3;
	allReflections[10][3] = 7;
	allReflections[10][7] = 4;
	allReflections[10][4] = 1;

	//S7-1-mmm
	allReflections[11][0] = 8;
	allReflections[11][8] = 0;
	allReflections[11][1] = 4;
	allReflections[11][4] = 7;
	allReflections[11][7] = 3;
	allReflections[11][3] = 6;
	allReflections[11][6] = 2;
	allReflections[11][2] = 5;
	allReflections[11][5] = 1;

	//S7-1-mm
	allReflections[12][0] = 8;
	allReflections[12][8] = 0;
	allReflections[12][1] = 3;
	allReflections[12][3] = 5;
	allReflections[12][5] = 7;
	allReflections[12][7] = 2;
	allReflections[12][2] = 4;
	allReflections[12][4] = 6;
	allReflections[12][6] = 1;

	//S7-1-m
	allReflections[13][0] = 8;
	allReflections[13][8] = 0;
	allReflections[13][1] = 2;
	allReflections[13][2] = 3;
	allReflections[13][3] = 4;
	allReflections[13][4] = 5;
	allReflections[13][5] = 6;
	allReflections[13][6] = 7;
	allReflections[13][7] = 1;

}


string Geometries::geometry9HBPYSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C7-1-p";
			break;
		case 1:
			return "C7-1-pp";
			break;
		case 2:
			return "C7-1-ppp";
			break;
		case 3:
			return "C7-1-mmm";
			break;
		case 4:
			return "C7-1-mm";
			break;
		case 5:
			return "C7-1-m";
			break;
		case 6:
			return "C2-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;
		case 7:
			return "C2-2 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;
		case 8:
			return "C2-3 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;
		case 9:
			return "C2-4 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;
		case 10:
			return "C2-5 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;
		case 11:
			return "C2-6 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;
		case 12:
			return "C2-7 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C7-1-p C7-1-pp C7-1-ppp C7-1-m C7-1-mm C7-1-mmm ) ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";
		case 3:
			return "P-4 ";
		case 4:
			return "P-5 ";
		case 5:
			return "P-6 ";
		case 6:
			return "P-7 ";
		case 7:
			return "P-8 ";
		case 8:
			return "S7-1-p";
		case 9:
			return "S7-1-pp";
		case 10:
			return "S7-1-ppp";
		case 11:
			return "S7-1-mmm";
		case 12:
			return "S7-1-mm";
		case 13:
			return "S7-1-m";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}



void Geometries::geometry9JTCotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(3);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][4] = 5;
	allReflections[0][5] = 4;
	allReflections[0][0] = 3;
	allReflections[0][3] = 0;
	allReflections[0][6] = 7;
	allReflections[0][7] = 6;
	allReflections[0][1] = 2;
	allReflections[0][2] = 1;

	//P-2
	allReflections[1][0] = 1;
	allReflections[1][1] = 0;
	allReflections[1][2] = 5;
	allReflections[1][5] = 2;
	allReflections[1][7] = 8;
	allReflections[1][8] = 7;
	allReflections[1][4] = 3;
	allReflections[1][3] = 4;

	//P-3
	allReflections[2][2] = 3;
	allReflections[2][3] = 2;
	allReflections[2][1] = 4;
	allReflections[2][4] = 1;
	allReflections[2][6] = 8;
	allReflections[2][8] = 6;
	allReflections[2][0] = 5;
	allReflections[2][5] = 0;

}


string Geometries::geometry9JTCSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
			break;
		case 1:
			return "C3-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}


void Geometries::geometry9JTDICotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 9;
	allReflections.resize(3);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}
	//P-1
	allReflections[0][1] = 7;
	allReflections[0][7] = 1;
	allReflections[0][4] = 8;
	allReflections[0][8] = 4;
	allReflections[0][3] = 6;
	allReflections[0][6] = 3;

	//P-2
	allReflections[1][1] = 2;
	allReflections[1][2] = 1;
	allReflections[1][0] = 3;
	allReflections[1][3] = 0;
	allReflections[1][4] = 5;
	allReflections[1][5] = 4;

	//P-3
	allReflections[2][2] = 7;
	allReflections[2][7] = 2;
	allReflections[2][0] = 6;
	allReflections[2][6] = 0;
	allReflections[2][5] = 8;
	allReflections[2][8] = 5;

}

string Geometries::geometry9JTDICSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C3-1-p";
			break;
		case 1:
			return "C3-1-m";
			break;

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ";
		case 1:
			return "P-2 ";
		case 2:
			return "P-3 ";

		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	return "ERROR";
}





///////////////////////////////////////////////////////////////////////////////////////////////////
void Geometries::geometry6OCotherSymmetries(
	std::vector< std::vector<int> > &allReflections)
{
	int size = 6;
	allReflections.resize(24);
	for (size_t i = 0; i < allReflections.size(); i++)
	{
		allReflections[i].resize(size);
		for (size_t j = 0; j < allReflections[i].size(); j++)
		{
			allReflections[i][j] = j;
		}
	}

	//P-1 ( C4-1-p C4-1-m C2-1 ) [h-0-5]
	allReflections[0][0] = 5;
	allReflections[0][5] = 0;

	//P-2  ( C2-8 )  [v(1)-0-5]
	allReflections[1][1] = 4;
	allReflections[1][4] = 1;
	allReflections[1][2] = 3;
	allReflections[1][3] = 2;

	//P-3  ( C2-9 )  [v(2)-0-5]
	allReflections[2][1] = 2;
	allReflections[2][2] = 1;
	allReflections[2][3] = 4;
	allReflections[2][4] = 3;

	//P-4  ( C4-2-p C4-2-m C2-2 )  [h-1-3]
	allReflections[3][1] = 3;
	allReflections[3][3] = 1;

	//P-5  ( C2-5 )  [v(1)-1-3]
	allReflections[4][2] = 5;
	allReflections[4][5] = 2;
	allReflections[4][0] = 4;
	allReflections[4][4] = 0;

	//P-6  ( C2-7 )  [v(2)-1-3]
	allReflections[5][0] = 2;
	allReflections[5][2] = 0;
	allReflections[5][4] = 5;
	allReflections[5][5] = 4;

	//P-7  ( C4-3-p C4-3-m C2-3 )  [h-2-4]
	allReflections[6][2] = 4;
	allReflections[6][4] = 2;

	//P-8  ( C2-6 )  [v(1)-2-4]
	allReflections[7][0] = 1;
	allReflections[7][1] = 0;
	allReflections[7][3] = 5;
	allReflections[7][5] = 3;

	//P-9  ( C2-4 )  [v(2)-2-4]
	allReflections[8][0] = 3;
	allReflections[8][3] = 0;
	allReflections[8][1] = 5;
	allReflections[8][5] = 1;

	//Inv
	allReflections[9][0] = 5;
	allReflections[9][5] = 0;
	allReflections[9][1] = 3;
	allReflections[9][3] = 1;
	allReflections[9][2] = 4;
	allReflections[9][4] = 2;

	//S6-1-p  [0 1 2]
	allReflections[10][0] = 3;
	allReflections[10][3] = 2;
	allReflections[10][2] = 5;
	allReflections[10][5] = 1;
	allReflections[10][1] = 4;
	allReflections[10][4] = 0;
	//S6-1-m  [0 1 2]
	allReflections[11][0] = 4;
	allReflections[11][4] = 1;
	allReflections[11][1] = 5;
	allReflections[11][5] = 2;
	allReflections[11][2] = 3;
	allReflections[11][3] = 1;

	//S6-2-p  [0 3 4]
	allReflections[12][0] = 1;
	allReflections[12][1] = 4;
	allReflections[12][4] = 5;
	allReflections[12][5] = 3;
	allReflections[12][3] = 2;
	allReflections[12][2] = 0;
	//S6-2-m  [0 3 4]
	allReflections[13][0] = 2;
	allReflections[13][2] = 3;
	allReflections[13][3] = 5;
	allReflections[13][5] = 4;
	allReflections[13][4] = 1;
	allReflections[13][1] = 0;

	//S6-3-p  [0 2 3]
	allReflections[14][0] = 4;
	allReflections[14][4] = 3;
	allReflections[14][3] = 5;
	allReflections[14][5] = 2;
	allReflections[14][2] = 1;
	allReflections[14][1] = 4;
	//S6-3-m  [0 2 3]
	allReflections[15][0] = 1;
	allReflections[15][1] = 2;
	allReflections[15][2] = 5;
	allReflections[15][5] = 3;
	allReflections[15][3] = 4;
	allReflections[15][4] = 0;

	//S6-4-p  [0 1 4]
	allReflections[16][0] = 2;
	allReflections[16][2] = 1;
	allReflections[16][1] = 5;
	allReflections[16][5] = 4;
	allReflections[16][4] = 3;
	allReflections[16][3] = 0;
	//S6-4-m  [0 1 4]
	allReflections[17][0] = 3;
	allReflections[17][3] = 4;
	allReflections[17][4] = 5;
	allReflections[17][5] = 1;
	allReflections[17][1] = 2;
	allReflections[17][2] = 0;

	//S4-1-p  [0 5]
	allReflections[18][0] = 5;
	allReflections[18][5] = 0;
	allReflections[18][2] = 1;
	allReflections[18][1] = 4;
	allReflections[18][4] = 3;
	allReflections[18][3] = 2;
	//S4-1-m  [0 5]
	allReflections[19][0] = 5;
	allReflections[19][5] = 0;
	allReflections[19][2] = 3;
	allReflections[19][3] = 4;
	allReflections[19][4] = 1;
	allReflections[19][1] = 2;

	//S4-2-p  [1 3]
	allReflections[20][3] = 1;
	allReflections[20][1] = 3;
	allReflections[20][0] = 4;
	allReflections[20][4] = 5;
	allReflections[20][5] = 2;
	allReflections[20][2] = 0;
	//S4-2-m  [1 3]
	allReflections[21][3] = 1;
	allReflections[21][1] = 3;
	allReflections[21][0] = 2;
	allReflections[21][2] = 5;
	allReflections[21][5] = 4;
	allReflections[21][4] = 0;

	//S4-3-p  [2 4]
	allReflections[22][2] = 4;
	allReflections[22][4] = 2;
	allReflections[22][0] = 1;
	allReflections[22][1] = 5;
	allReflections[22][5] = 3;
	allReflections[22][3] = 0;
	//S4-3-m  [2 4]
	allReflections[23][2] = 4;
	allReflections[23][4] = 2;
	allReflections[23][0] = 3;
	allReflections[23][3] = 5;
	allReflections[23][5] = 1;
	allReflections[23][1] = 0;

}

string Geometries::geometry6OCSymmetryFlags(
	int iSymmetry,
	int symmetryType)
{
	if (symmetryType == 0)
	{
		switch (iSymmetry)
		{
		case 0:
			return "C4-1-p";
		case 1:
			return "C2-1 ( C4-2-p C4-2-m C2-2 C4-3-p C4-3-m C2-3 C2-8 C2-9 )";
		case 2:
			return "C4-1-m";
		case 3:
			return "C4-2-p";
		case 4:
			return "C2-2 ( C4-1-p C4-1-m C2-1 C4-3-p C4-3-m C2-3 C2-5 C2-7 )";
		case 5:
			return "C4-2-m";
		case 6:
			return "C4-3-p";
		case 7:
			return "C2-3 ( C4-1-p C4-1-m C2-1 C4-2-p C4-2-m C2-2 C2-4 C2-6 )";
		case 8:
			return "C4-3-m";
		case 9:
			return "C2-4 ( C4-3-p C4-3-m C2-3 C3-2-p C3-2-m C3-3-p C3-3-m C2-6 )";
		case 10:
			return "C2-5 ( C4-2-p C4-2-m C2-2 C3-3-p C3-3-m C3-4-p C3-4-m C2-7 )";
		case 11:
			return "C2-6 ( C4-3-p C4-3-m C2-3 C3-1-p C3-1-m C3-4-p C3-4-m C2-4 )";
		case 12:
			return "C2-7  ( C4-2-p C4-2-m C2-2 C3-1-p C3-1-m C3-2-p C3-2-m C2-5 )";
		case 13:
			return "C2-8  ( C4-1-p C4-1-m C2-1 C3-2-p C3-2-m C3-4-p C3-4-m C2-9 )";
		case 14:
			return "C2-9  ( C4-1-p C4-1-m C2-1 C3-1-p C3-1-m C3-3-p C3-3-m C2-8 )";
		case 15:
			return "C3-1-p";
		case 16:
			return "C3-1-m";
		case 17:
			return "C3-2-p";
		case 18:
			return "C3-2-m";
		case 19:
			return "C3-3-p";
		case 20:
			return "C3-3-m";
		case 21:
			return "C3-4-p";
		case 22:
			return "C3-4-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}
	}
	else
	{
		switch (iSymmetry)
		{
		case 0:
			return "P-1 ( C4-1-p C4-1-m C2-1 )";
		case 1:
			return "P-2  ( C2-8 )";
		case 2:
			return "P-3  ( C2-9 )";
		case 3:
			return "P-4  ( C4-2-p C4-2-m C2-2 )";
		case 4:
			return "P-5  ( C2-5 )";
		case 5:
			return "P-6  ( C2-7 )";
		case 6:
			return "P-7  ( C4-3-p C4-3-m C2-3 )";
		case 7:
			return "P-8  ( C2-6 )";
		case 8:
			return "P-9  ( C2-4 )";
		case 9:
			return "Inv";
		case 10:
			return "S6-1-p";
		case 11:
			return "S6-1-m";
		case 12:
			return "S6-2-p";
		case 13:
			return "S6-2-m";
		case 14:
			return "S6-3-p";
		case 15:
			return "S6-3-m";
		case 16:
			return "S6-4-p";
		case 17:
			return "S6-4-m";
		case 18:
			return "S4-1-p";
		case 19:
			return "S4-1-m";
		case 20:
			return "S4-2-p";
		case 21:
			return "S4-2-m";
		case 22:
			return "S4-3-p";
		case 23:
			return "S4-3-m";
		default:
			cout << "rotation on CauchyIndex::rotationString not found" << endl;
			exit(1);
			break;
		}

	}

	return "ERROR";
}


