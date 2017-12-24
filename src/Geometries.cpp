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

	case 50:
		return geometry5TBPY(mol0, cutAngle, reflectionOperation);
		break;

	case 51:
		return geometry5SPY(mol0, cutAngle, reflectionOperation);
		break;

	case 52:
		return geometry5VOC(mol0, cutAngle, reflectionOperation);
		break;

	case 60:
		return geometry6OC(mol0, cutAngle, reflectionOperation);
		break;

	case 61:
		return geometry6TPR(mol0, cutAngle, reflectionOperation);
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

	case 90:
		return geometry9TCTPR(mol0, cutAngle, reflectionOperation);
		break;

	case 91:
		return geometry9CSAPR(mol0, cutAngle, reflectionOperation);
		break;

	case 92:
		return geometry9MFF(mol0, cutAngle, reflectionOperation);
		break;

	case 100:
		return geometry10JMBIC(mol0, cutAngle, reflectionOperation);
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

std::vector<double> Geometries::geometry4Tetrahedron(
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
	//c2-0
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi;
	// c2-1
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	// c2-2
	vectorRotations[8] = mol0[1].x;
	vectorRotations[9] = mol0[1].y;
	vectorRotations[10] = mol0[1].z;
	vectorRotations[11] = auxMath_._pi;
	// c2-3
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;
	// c4-1
	vectorRotations[16] = 0.00000000000000000;
	vectorRotations[17] = 0.00000000000000000;
	vectorRotations[18] = 1.00000000000000000;
	vectorRotations[19] = auxMath_._pi / 2.0e0;
	// c4-2
	vectorRotations[20] = 0.00000000000000000;
	vectorRotations[21] = 0.00000000000000000;
	vectorRotations[22] = 1.00000000000000000;
	vectorRotations[23] = auxMath_._pi;
	// c4-3
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
	//c2-0
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi;
	// c2-1
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	// c2-2
	vectorRotations[8] = mol0[1].x;
	vectorRotations[9] = mol0[1].y;
	vectorRotations[10] = mol0[1].z;
	vectorRotations[11] = auxMath_._pi;
	// c2-3
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;
	// c4-1
	vectorRotations[16] = 0.00000000000000000;
	vectorRotations[17] = 0.00000000000000000;
	vectorRotations[18] = 1.00000000000000000;
	vectorRotations[19] = auxMath_._pi / 2.0e0;
	// c4-2
	vectorRotations[20] = 0.00000000000000000;
	vectorRotations[21] = 0.00000000000000000;
	vectorRotations[22] = 1.00000000000000000;
	vectorRotations[23] = auxMath_._pi;
	// c4-3
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


std::vector<double> Geometries::geometry5TBPY(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 5;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = 0.9 * auxMath_._pi;
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

	//c3 - 1
	vectorRotations[0] = 1.0e0;
	vectorRotations[1] = 0.0e0;
	vectorRotations[2] = 0.0e0;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 2
	vectorRotations[4] = 1.0e0;
	vectorRotations[5] = 0.0e0;
	vectorRotations[6] = 0.0e0;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;
	//c2 - 1
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = auxMath_._pi;
	//c2 - 2
	vectorRotations[12] = mol0[3].x;
	vectorRotations[13] = mol0[3].y;
	vectorRotations[14] = mol0[3].z;
	vectorRotations[15] = auxMath_._pi;
	//c2 - 3
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
	cutAngle = 0.75 * auxMath_._pi;
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

	// 1 c3 por face.
	// tem 3 pontos que geram os c4 e c3.
	// 6 c2 cruzam as 12 arestas.
	//c4 - 1 (ref 0)
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//c4 - 2 (ref 0)
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = auxMath_._pi;
	//c4 - 3 (ref 0)
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 3.0e0 * auxMath_._pi / 2.0e0;

	//c4 - 1 (ref 1)
	vectorRotations[12] = mol0[1].x;
	vectorRotations[13] = mol0[1].y;
	vectorRotations[14] = mol0[1].z;
	vectorRotations[15] = auxMath_._pi / 2.0e0;
	//c4 - 2 (ref 1)
	vectorRotations[16] = mol0[1].x;
	vectorRotations[17] = mol0[1].y;
	vectorRotations[18] = mol0[1].z;
	vectorRotations[19] = auxMath_._pi;
	//c4 - 3 (ref 1)
	vectorRotations[20] = mol0[1].x;
	vectorRotations[21] = mol0[1].y;
	vectorRotations[22] = mol0[1].z;
	vectorRotations[23] = 3.0e0 * auxMath_._pi / 2.0e0;

	//c4 - 1 (ref 2)
	vectorRotations[24] = mol0[2].x;
	vectorRotations[25] = mol0[2].y;
	vectorRotations[26] = mol0[2].z;
	vectorRotations[27] = auxMath_._pi / 2.0e0;
	//c4 - 2 (ref 2)
	vectorRotations[28] = mol0[2].x;
	vectorRotations[29] = mol0[2].y;
	vectorRotations[30] = mol0[2].z;
	vectorRotations[31] = auxMath_._pi;
	//c4 - 3 (ref 2)
	vectorRotations[32] = mol0[2].x;
	vectorRotations[33] = mol0[2].y;
	vectorRotations[34] = mol0[2].z;
	vectorRotations[35] = 3.0e0 * auxMath_._pi / 2.0e0;

	//c2 - 0 - 1
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[36] = auxReferenceAxis[0];
	vectorRotations[37] = auxReferenceAxis[1];
	vectorRotations[38] = auxReferenceAxis[2];
	vectorRotations[39] = auxMath_._pi;

	//c2 - 0 - 2
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[40] = auxReferenceAxis[0];
	vectorRotations[41] = auxReferenceAxis[1];
	vectorRotations[42] = auxReferenceAxis[2];
	vectorRotations[43] = auxMath_._pi;

	//c2 - 0 - 3
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[44] = auxReferenceAxis[0];
	vectorRotations[45] = auxReferenceAxis[1];
	vectorRotations[46] = auxReferenceAxis[2];
	vectorRotations[47] = auxMath_._pi;

	//c2 - 0 - 4
	auxReferenceAxis[0] = 0.5e0 * (mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[48] = auxReferenceAxis[0];
	vectorRotations[49] = auxReferenceAxis[1];
	vectorRotations[50] = auxReferenceAxis[2];
	vectorRotations[51] = auxMath_._pi;

	//c2 - 1 - 2
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[52] = auxReferenceAxis[0];
	vectorRotations[53] = auxReferenceAxis[1];
	vectorRotations[54] = auxReferenceAxis[2];
	vectorRotations[55] = auxMath_._pi;

	//c2 - 1 - 4
	auxReferenceAxis[0] = 0.5e0 * (mol0[1].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[1].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[1].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[56] = auxReferenceAxis[0];
	vectorRotations[57] = auxReferenceAxis[1];
	vectorRotations[58] = auxReferenceAxis[2];
	vectorRotations[59] = auxMath_._pi;

	//c3 ---- 0 - 1 - 2
	auxReferenceAxis[0] = (mol0[0].x + mol0[1].x + mol0[2].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[1].y + mol0[2].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[1].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[60] = auxReferenceAxis[0];
	vectorRotations[61] = auxReferenceAxis[1];
	vectorRotations[62] = auxReferenceAxis[2];
	vectorRotations[63] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[64] = auxReferenceAxis[0];
	vectorRotations[65] = auxReferenceAxis[1];
	vectorRotations[66] = auxReferenceAxis[2];
	vectorRotations[67] = 4.0e0 * auxMath_._pi / 3.0e0;

	//c3 ---- 0 - 2 - 3
	auxReferenceAxis[0] = (mol0[0].x + mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[68] = auxReferenceAxis[0];
	vectorRotations[69] = auxReferenceAxis[1];
	vectorRotations[70] = auxReferenceAxis[2];
	vectorRotations[71] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[72] = auxReferenceAxis[0];
	vectorRotations[73] = auxReferenceAxis[1];
	vectorRotations[74] = auxReferenceAxis[2];
	vectorRotations[75] = 4.0e0 * auxMath_._pi / 3.0e0;

	//c3 ---- 0 - 3 - 4
	auxReferenceAxis[0] = (mol0[0].x + mol0[3].x + mol0[4].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[3].y + mol0[4].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[3].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[76] = auxReferenceAxis[0];
	vectorRotations[77] = auxReferenceAxis[1];
	vectorRotations[78] = auxReferenceAxis[2];
	vectorRotations[79] = 2.0e0 * auxMath_._pi / 3.0e0;
	vectorRotations[80] = auxReferenceAxis[0];
	vectorRotations[81] = auxReferenceAxis[1];
	vectorRotations[82] = auxReferenceAxis[2];
	vectorRotations[83] = 4.0e0 * auxMath_._pi / 3.0e0;

	//c3 ---- 0 - 4 - 1
	auxReferenceAxis[0] = (mol0[0].x + mol0[4].x + mol0[1].x);
	auxReferenceAxis[1] = (mol0[0].y + mol0[4].y + mol0[1].y);
	auxReferenceAxis[2] = (mol0[0].z + mol0[4].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[84] = auxReferenceAxis[0];
	vectorRotations[85] = auxReferenceAxis[1];
	vectorRotations[86] = auxReferenceAxis[2];
	vectorRotations[87] = 2.0e0 * auxMath_._pi / 3.0e0;
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

	// 1 c3 por face.
	// tem 3 pontos que geram os c4 e c3.
	// 6 c2 cruzam as 12 arestas.
	//c3 - 1 (ref 0)
	vectorRotations[0] = 0.000000000;
	vectorRotations[1] = 0.000000000;
	vectorRotations[2] = 1.000000000;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 2 (ref 0)
	vectorRotations[4] = 0.000000000;
	vectorRotations[5] = 0.000000000;
	vectorRotations[6] = 1.000000000;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 3.0e0;

	//c2-1
	auxReferenceAxis[0] = 0.5e0 * (mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = auxMath_._pi;

	//c2-2
	auxReferenceAxis[0] = 0.5e0 * (mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 * (mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 * (mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi;

	//c2-2
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


std::vector<double> Geometries::geometry7COC(
	std::vector<CoordXYZ> &mol0,
	double & cutAngle,
	std::vector<int> &reflectionOperation)
{
	int size = 7;
	mol0.resize(size);
	reflectionOperation.resize(size);
	cutAngle = auxMath_._pi / 2.0e0;;
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

	//c3 - 1
	vectorRotations[0] = mol0[5].x;
	vectorRotations[1] = mol0[5].y;
	vectorRotations[2] = mol0[5].z;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 2
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

	//c5 - 1
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = 2.0e0 * auxMath_._pi / 5.0e0;
	//c5 - 2
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 4.0e0 * auxMath_._pi / 5.0e0;
	//c5 - 3
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = 6.0e0 * auxMath_._pi / 5.0e0;
	//c5 - 4
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = 8.0e0 * auxMath_._pi / 5.0e0;

	//c2 - 1
	vectorRotations[16] = mol0[2].x;
	vectorRotations[17] = mol0[2].y;
	vectorRotations[18] = mol0[2].z;
	vectorRotations[19] = auxMath_._pi;
	//c2 - 2
	vectorRotations[20] = mol0[3].x;
	vectorRotations[21] = mol0[3].y;
	vectorRotations[22] = mol0[3].z;
	vectorRotations[23] = auxMath_._pi;
	//c2 - 3
	vectorRotations[24] = mol0[4].x;
	vectorRotations[25] = mol0[4].y;
	vectorRotations[26] = mol0[4].z;
	vectorRotations[27] = auxMath_._pi;
	//c2 - 4
	vectorRotations[28] = mol0[5].x;
	vectorRotations[29] = mol0[5].y;
	vectorRotations[30] = mol0[5].z;
	vectorRotations[31] = auxMath_._pi;
	//c2 - 5
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

	//c2-1
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
	//c4 - 1
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//c4 - 2 (c2)
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	//c4 - 3
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = 3.0e0 * auxMath_._pi / 2.0e0;
	// c2 - 1 e 7
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
	// c2 - 7 e 5
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
	// c2 - 5 e 6
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
	// c2 - 6 e 2
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
	cutAngle = 2.0e0 * auxMath_._pi / 3.0e0;;
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

	// c2 -z
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi;
	// c2 -1
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	// c2 -2
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

	// c2 -z
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

	//c6 - 1
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 3.0e0;
	//c6 - 2
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c6 - 3
	vectorRotations[8] = mol0[0].x;
	vectorRotations[9] = mol0[0].y;
	vectorRotations[10] = mol0[0].z;
	vectorRotations[11] = auxMath_._pi;
	//c6 - 4
	vectorRotations[12] = mol0[0].x;
	vectorRotations[13] = mol0[0].y;
	vectorRotations[14] = mol0[0].z;
	vectorRotations[15] = -2.0e0 * auxMath_._pi / 3.0e0;
	//c6 - 5
	vectorRotations[16] = mol0[0].x;
	vectorRotations[17] = mol0[0].y;
	vectorRotations[18] = mol0[0].z;
	vectorRotations[19] = -auxMath_._pi / 3.0e0;
	//c2 - 1
	vectorRotations[20] = mol0[2].x;
	vectorRotations[21] = mol0[2].y;
	vectorRotations[22] = mol0[2].z;
	vectorRotations[23] = auxMath_._pi;
	//c2 - 2
	vectorRotations[24] = mol0[3].x;
	vectorRotations[25] = mol0[3].y;
	vectorRotations[26] = mol0[3].z;
	vectorRotations[27] = auxMath_._pi;
	//c2 - 3
	vectorRotations[28] = mol0[4].x;
	vectorRotations[29] = mol0[4].y;
	vectorRotations[30] = mol0[4].z;
	vectorRotations[31] = auxMath_._pi;
	//c2l - 1
	auxReferenceAxis[0] = 0.5e0 *(mol0[4].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[4].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[4].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = auxMath_._pi;
	//c2l - 2
	auxReferenceAxis[0] = 0.5e0 *(mol0[2].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[2].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[2].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[36] = auxReferenceAxis[0];
	vectorRotations[37] = auxReferenceAxis[1];
	vectorRotations[38] = auxReferenceAxis[2];
	vectorRotations[39] = auxMath_._pi;
	//c2l - 3
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

	// C4s
	// c4f1 - 0-6 - 1
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[0] = auxReferenceAxis[0];
	vectorRotations[1] = auxReferenceAxis[1];
	vectorRotations[2] = auxReferenceAxis[2];
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	// c4f1 - 0-6 - 2
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[4] = auxReferenceAxis[0];
	vectorRotations[5] = auxReferenceAxis[1];
	vectorRotations[6] = auxReferenceAxis[2];
	vectorRotations[7] = auxMath_._pi;
	// c4f1 - 0-6 - 3
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[6].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[6].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[8] = auxReferenceAxis[0];
	vectorRotations[9] = auxReferenceAxis[1];
	vectorRotations[10] = auxReferenceAxis[2];
	vectorRotations[11] = -auxMath_._pi / 2.0e0;
	// c4f2 - 0-5 - 1
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = auxMath_._pi / 2.0e0;
	// c4f2 - 0-5 - 2
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[16] = auxReferenceAxis[0];
	vectorRotations[17] = auxReferenceAxis[1];
	vectorRotations[18] = auxReferenceAxis[2];
	vectorRotations[19] = auxMath_._pi;
	// c4f2 - 0-5 - 3
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[5].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[5].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[5].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[20] = auxReferenceAxis[0];
	vectorRotations[21] = auxReferenceAxis[1];
	vectorRotations[22] = auxReferenceAxis[2];
	vectorRotations[23] = -auxMath_._pi / 2.0e0;
	// c4f3 - 0-4 - 1
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[24] = auxReferenceAxis[0];
	vectorRotations[25] = auxReferenceAxis[1];
	vectorRotations[26] = auxReferenceAxis[2];
	vectorRotations[27] = auxMath_._pi / 2.0e0;
	// c4f3 - 0-4 - 2
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[28] = auxReferenceAxis[0];
	vectorRotations[29] = auxReferenceAxis[1];
	vectorRotations[30] = auxReferenceAxis[2];
	vectorRotations[31] = auxMath_._pi;
	// c4f3 - 0-4 - 3
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[32] = auxReferenceAxis[0];
	vectorRotations[33] = auxReferenceAxis[1];
	vectorRotations[34] = auxReferenceAxis[2];
	vectorRotations[35] = -auxMath_._pi / 2.0e0;
	//c3 - 0 - 1
	vectorRotations[36] = mol0[0].x;
	vectorRotations[37] = mol0[0].y;
	vectorRotations[38] = mol0[0].z;
	vectorRotations[39] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 0 - 2
	vectorRotations[40] = mol0[0].x;
	vectorRotations[41] = mol0[0].y;
	vectorRotations[42] = mol0[0].z;
	vectorRotations[43] = 4.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 1 - 1
	vectorRotations[44] = mol0[1].x;
	vectorRotations[45] = mol0[1].y;
	vectorRotations[46] = mol0[1].z;
	vectorRotations[47] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 1 - 2
	vectorRotations[48] = mol0[1].x;
	vectorRotations[49] = mol0[1].y;
	vectorRotations[50] = mol0[1].z;
	vectorRotations[51] = 4.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 2 - 1
	vectorRotations[52] = mol0[2].x;
	vectorRotations[53] = mol0[2].y;
	vectorRotations[54] = mol0[2].z;
	vectorRotations[55] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 2 - 2
	vectorRotations[56] = mol0[2].x;
	vectorRotations[57] = mol0[2].y;
	vectorRotations[58] = mol0[2].z;
	vectorRotations[59] = 4.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 4 - 1
	vectorRotations[60] = mol0[4].x;
	vectorRotations[61] = mol0[4].y;
	vectorRotations[62] = mol0[4].z;
	vectorRotations[63] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 4 - 2
	vectorRotations[64] = mol0[4].x;
	vectorRotations[65] = mol0[4].y;
	vectorRotations[66] = mol0[4].z;
	vectorRotations[67] = 4.0e0 * auxMath_._pi / 3.0e0;
	// c2 - aresta - 1
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[1].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[1].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[1].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[68] = auxReferenceAxis[0];
	vectorRotations[69] = auxReferenceAxis[1];
	vectorRotations[70] = auxReferenceAxis[2];
	vectorRotations[71] = auxMath_._pi;
	// c2 - aresta - 2
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[2].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[2].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[2].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[72] = auxReferenceAxis[0];
	vectorRotations[73] = auxReferenceAxis[1];
	vectorRotations[74] = auxReferenceAxis[2];
	vectorRotations[75] = auxMath_._pi;
	// c2 - aresta - 3
	auxReferenceAxis[0] = 0.5e0 *(mol0[0].x + mol0[3].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[0].y + mol0[3].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[0].z + mol0[3].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[76] = auxReferenceAxis[0];
	vectorRotations[77] = auxReferenceAxis[1];
	vectorRotations[78] = auxReferenceAxis[2];
	vectorRotations[79] = auxMath_._pi;
	// c2 - aresta - 4
	auxReferenceAxis[0] = 0.5e0 *(mol0[2].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[2].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[2].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[80] = auxReferenceAxis[0];
	vectorRotations[81] = auxReferenceAxis[1];
	vectorRotations[82] = auxReferenceAxis[2];
	vectorRotations[83] = auxMath_._pi;
	// c2 - aresta - 5
	auxReferenceAxis[0] = 0.5e0 *(mol0[1].x + mol0[4].x);
	auxReferenceAxis[1] = 0.5e0 *(mol0[1].y + mol0[4].y);
	auxReferenceAxis[2] = 0.5e0 *(mol0[1].z + mol0[4].z);
	auxMath_.normalize(auxReferenceAxis);
	vectorRotations[84] = auxReferenceAxis[0];
	vectorRotations[85] = auxReferenceAxis[1];
	vectorRotations[86] = auxReferenceAxis[2];
	vectorRotations[87] = auxMath_._pi;
	// c2 - aresta - 6
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

	//c2 - 1
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi;
	//c2 - 2
	vectorRotations[4] = mol0[6].x;
	vectorRotations[5] = mol0[6].y;
	vectorRotations[6] = mol0[6].z;
	vectorRotations[7] = auxMath_._pi;
	//c2 - 3
	vectorRotations[8] = mol0[5].x;
	vectorRotations[9] = mol0[5].y;
	vectorRotations[10] = mol0[5].z;
	vectorRotations[11] = auxMath_._pi;
	// produto vetorial
	auxReferenceAxis = auxMath_.vectorProduct(
		mol0[0].x, mol0[0].y, mol0[0].z,
		mol0[6].x, mol0[6].y, mol0[6].z);
	auxMath_.normalize(auxReferenceAxis);
	//c3 - 1
	vectorRotations[12] = auxReferenceAxis[0];
	vectorRotations[13] = auxReferenceAxis[1];
	vectorRotations[14] = auxReferenceAxis[2];
	vectorRotations[15] = 2.0e0 * auxMath_._pi / 3.0e0;
	//c3 - 2
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
	vector<double> auxReferenceAxis(3);

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

	//c4 - 1
	vectorRotations[0] = mol0[0].x;
	vectorRotations[1] = mol0[0].y;
	vectorRotations[2] = mol0[0].z;
	vectorRotations[3] = auxMath_._pi / 2.0e0;
	//c4 - 2
	vectorRotations[4] = mol0[0].x;
	vectorRotations[5] = mol0[0].y;
	vectorRotations[6] = mol0[0].z;
	vectorRotations[7] = auxMath_._pi;
	//c4 - 3
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
	vector<double> auxReferenceAxis(3);

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

	return vectorRotations;
}






std::vector<double> Geometries::geometry10JMBIC(
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








