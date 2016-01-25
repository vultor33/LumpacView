#include "Ligand.h"

#include <fstream>

#include "Coordstructs.h"
#include "AuxMath.h"

using namespace std;

Ligand::Ligand() {}

Ligand::~Ligand() {}

// xm1 xm2 - x1 e x2 monodentado
// xb1 xb2 xb3 - 1 e 2 sao os atomos, 3 e a direcao.
// xt1 xt2 xt3 xt4 - 1 2 e 3 sao os atomos. 4 e pra saber a direcao que ele nao pode ir.
bool Ligand::initializeLigand(vector<CoordXYZ> &coord_in)
{
	int size = coord_in.size();
	if (
		(coord_in[size - 1].atomlabel == "Xm2") &&
		(coord_in[size - 2].atomlabel == "Xm1"))
		chelation = 1;
	else if (
		(coord_in[size - 1].atomlabel == "Xb3") &&
		(coord_in[size - 2].atomlabel == "Xb2") &&
		(coord_in[size - 3].atomlabel == "Xb1"))
		chelation = 2;
	else if (
		(coord_in[size - 1].atomlabel == "Xt4") &&
		(coord_in[size - 2].atomlabel == "Xt3") &&
		(coord_in[size - 3].atomlabel == "Xt2") &&
		(coord_in[size - 4].atomlabel == "Xt1"))
		chelation = 3;
	else
		return false;

	CoordXYZ X1;
	X1.atomlabel = "X1";
	CoordXYZ X2;
	X2.atomlabel = "X2";
	if (chelation == 1)
	{
		X1.x = coord_in[size - 2].x;
		X1.y = coord_in[size - 2].y;
		X1.z = coord_in[size - 2].z;

		X2.x = coord_in[size - 1].x;
		X2.y = coord_in[size - 1].y;
		X2.z = coord_in[size - 1].z;
	}
	else if (chelation == 2)
	{
		X1.x = 0.5e0 * (coord_in[size - 3].x + coord_in[size - 2].x);
		X1.y = 0.5e0 * (coord_in[size - 3].y + coord_in[size - 2].y);
		X1.z = 0.5e0 * (coord_in[size - 3].z + coord_in[size - 2].z);

		X2.x = coord_in[size - 1].x;
		X2.y = coord_in[size - 1].y;
		X2.z = coord_in[size - 1].z;
	}
	else if (chelation == 3)
	{
		AuxMath math_;
		vector<double> centroid = math_.triangleCentroid
			(coord_in[size - 4].x, coord_in[size - 4].y, coord_in[size - 4].z,
				coord_in[size - 3].x, coord_in[size - 3].y, coord_in[size - 3].z,
				coord_in[size - 2].x, coord_in[size - 2].y, coord_in[size - 2].z);
		X1.x = centroid[0];
		X1.y = centroid[1];
		X1.z = centroid[2];

		vector<double> normal = math_.normalVectorFrom3Points
			(coord_in[size - 4].x, coord_in[size - 4].y, coord_in[size - 4].z,
				coord_in[size - 3].x, coord_in[size - 3].y, coord_in[size - 3].z,
				coord_in[size - 2].x, coord_in[size - 2].y, coord_in[size - 2].z);

		X2.x = X1.x + normal[0];
		X2.y = X1.y + normal[1];
		X2.z = X1.z + normal[2];

		double escalar = math_.escalarProduct(
			X2.x, X2.y, X2.z,
			coord_in[size - 1].x, coord_in[size - 1].y, coord_in[size - 1].z);
		if (escalar < 0.0e0)
		{
			X2.x *= -1.0e0;
			X2.y *= -1.0e0;
			X2.z *= -1.0e0;
		}
	}
	else
		return false;

	size += 1 - chelation;
	coord.resize(size);
	for (int i = 0; i < size - 2; i++)
	{
		coord[i] = coord_in[i];

	}
	coord[size - 2] = X1;
	coord[size - 1] = X2;

#ifdef _DEBUG
	printXyzLigandDirection();
#endif

	return true;
}

void Ligand::printXyzLigandDirection()
{
	ofstream out_("xyzTeste.xyz");
	int size = coord.size();
	out_ << size << endl
		<< "useless" << endl;
	for (int i = 0; i < (size-2); i++)
	{
		out_ << coord[i].atomlabel << "   "
			<< coord[i].x << "   "
			<< coord[i].y << "   "
			<< coord[i].z << endl;
	}
	int x1 = size - 2;
	int x2 = size - 1;
	out_ << "Au  " << "   "
		<< coord[x1].x << "   "
		<< coord[x1].y << "   "
		<< coord[x1].z << endl;

	out_ << "Au " << "   "
		<< coord[x2].x << "   "
		<< coord[x2].y << "   "
		<< coord[x2].z << endl;

	out_.close();

}



