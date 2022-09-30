#include "MarquesEnantiomers.h"

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "Coordstructs.h"
#include "AuxMathGa.h"
#include "Hungarian.h"

using namespace std;
using namespace zerg;

MarquesEnantiomers::MarquesEnantiomers()
{
	tolRmsFit = 1.0e-4;
	tolIdentificalStructRmsd = 5.0e-2;
	nCycles = 100000;
	minRepeat = 50;
	elementDiff = 20.0e0;
	kCount = 0;
}

MarquesEnantiomers::~MarquesEnantiomers(){}

void MarquesEnantiomers::setSeed(int seed)
{
	AuxMathGa::set_seed(seed);
}

double MarquesEnantiomers::marquesRmsd(
	vector<CoordXYZ> &mol1,
	vector<CoordXYZ> &mol2)
{
	for (size_t i = 0; i < mol1.size(); i++)
	{
		mol1[i].mass = setMass(mol1[i].atomlabel);
		mol2[i].mass = setMass(mol2[i].atomlabel);
	}
	return assignStruct(mol1, mol2);
}


void MarquesEnantiomers::calculateMarquesEnantiomers(
	string fileXyz1, 
	string fileXyz2)
{
	vector<CoordXYZ> mol1 = readXyz(fileXyz1);
	vector<CoordXYZ> mol2 = readXyz(fileXyz2);
	vector<CoordXYZ> mol2Mirror = mol2;
	for (size_t i = 0; i < mol2.size(); i++)
	{
		mol2Mirror[i].y *= -1.0e0; //xz plane reflection
	}
	double rmsd1 = assignStruct(mol1, mol2);

	cout << "count:  " << kCount << endl;

	if (rmsd1 < tolIdentificalStructRmsd)
	{
		cout << "SAO IDENTICOS ESSA MERDA  " << rmsd1 << endl;
		return;
	}

	double rmsd2 = assignStruct(mol1, mol2Mirror);

	cout << "count:  " << kCount << endl;

	if ((rmsd2 < rmsd1) && (rmsd2 < tolIdentificalStructRmsd))
	{
		cout << "SAO ENANTIOMEROS" << endl;
	}
	else
		cout << "VIDA SEGUE" << endl;

}

double MarquesEnantiomers::assignStruct(
	vector<CoordXYZ> &mol1,
	vector<CoordXYZ> &mol2)
{
	translateToCenterOfMass(mol1);
	translateToCenterOfMass(mol2);
	vector<int> connect1(mol1.size());
	vector<int> connect2(mol2.size());
	vector<CoordXYZ> mol2Initial = mol2;

	int rmsd_repeat;
	double rmsdmin;
	vector<int> connectMin(mol1.size());
	vector<CoordXYZ> mol2Min;
	for (size_t i = 0; i < mol1.size(); i++)
	{
		connect1[i] = i;
		connect2[i] = i;
	}
	for (int k = 0; k < nCycles; k++)
	{
		kCount = k;
		
		mol2 = mol2Initial;

		eulerRotation(mol2);

		double rmsI = rmsFitConnect(mol1, connect1, mol2, connect2);

		vector< vector<double> > matrDist = distanceMatrix(mol1, mol2);

		HungarianAlgorithm assign_;

		double cost = assign_.Solve(matrDist, connect2);

		quatFit(mol1, connect1, mol2, connect2);

		double rmsd = rmsFitConnect(mol1, connect1, mol2, connect2);

		if (k == 0)
		{
			rmsd_repeat = 1;
			rmsdmin = rmsd;
			connectMin = connect2;
			mol2Min = mol2;
		}
		else
		{
			double rmsdTol = rmsd - rmsdmin;
			if (abs(rmsdTol) < tolRmsFit)
			{
				rmsd_repeat++;
				if (rmsd_repeat == minRepeat)
				{					
					kCount = k;
					return rmsdmin;
				}
			}
			else if (rmsdTol < 0.0e0)
			{
					rmsd_repeat = 1;
					rmsdmin = rmsd;
					connectMin = connect2;
					mol2Min = mol2;
			}
		}
	}
	cout << "WARNING - MARQUESRMSD DIDN'T CONVERGED" << endl;
	return rmsdmin;
}

vector<CoordXYZ> MarquesEnantiomers::readXyz(string xyzName)
{
	ifstream xyzFile_(xyzName.c_str());
	string line;
	stringstream convert;
	int natoms;
	getline(xyzFile_, line);
	convert << line;
	convert >> natoms;
	getline(xyzFile_, line);

	vector<CoordXYZ> mol(natoms);
	for (int i = 0; i < natoms; i++)
	{
		getline(xyzFile_, line);
		stringstream convert1;
		convert1 << line;
		convert1 >> mol[i].atomlabel
			>> mol[i].x
			>> mol[i].y
			>> mol[i].z;
		mol[i].mass = setMass(mol[i].atomlabel);
	}
	return mol;
}

void MarquesEnantiomers::translateToCenterOfMass(vector<CoordXYZ> &mol)
{
	double xcm = 0.0e0;
	double ycm = 0.0e0;
	double zcm = 0.0e0;
	double totalMass = 0.0e0;
	for (size_t i = 0; i < mol.size(); i++)
	{
		double massI = mol[i].mass;
		totalMass += massI;
		xcm += mol[i].x * massI;
		ycm += mol[i].y * massI;
		zcm += mol[i].z * massI;
	}
	xcm /= totalMass;
	ycm /= totalMass;
	zcm /= totalMass;
	for (size_t i = 0; i < mol.size(); i++)
	{
		mol[i].x -= xcm;
		mol[i].y -= ycm;
		mol[i].z -= zcm;
	}
}

void MarquesEnantiomers::eulerRotation(vector<CoordXYZ> &mol)
{
	double TWOPI = 2.0e0 * acos(-1.0e0);
	double upperLimit = 1.0e0;
	double RAND;
	RAND = AuxMathGa::randomNumber(0.0e0, upperLimit);
	double PHI = TWOPI*RAND;
	RAND = AuxMathGa::randomNumber(0.0e0, upperLimit);
	double CSTHTA = 2.0e0*RAND - 1.0e0;
	RAND = AuxMathGa::randomNumber(0.0e0, upperLimit);
	double CHI = TWOPI*RAND;
	double THTA = acos(CSTHTA);
	double SNTHTA = sin(THTA);
	double SNPHI = sin(PHI);
	double CSPHI = cos(PHI);
	double SNCHI = sin(CHI);
	double CSCHI = cos(CHI);
	double RXX = CSTHTA*CSPHI*CSCHI - SNPHI*SNCHI;
	double RXY = -CSTHTA*CSPHI*SNCHI - SNPHI*CSCHI;
	double RXZ = SNTHTA*CSPHI;
	double RYX = CSTHTA*SNPHI*CSCHI + CSPHI*SNCHI;
	double RYY = -CSTHTA*SNPHI*SNCHI + CSPHI*CSCHI;
	double RYZ = SNTHTA*SNPHI;
	double RZX = -SNTHTA*CSCHI;
	double RZY = SNTHTA*SNCHI;
	double RZZ = CSTHTA;
	vector<double> Q(3 * mol.size());
	int size = mol.size();
	for (size_t i = 0; i < mol.size(); i++)
	{
		Q[i] = mol[i].x*RXX + mol[i].y*RXY + mol[i].z*RXZ;
		Q[i + size] = mol[i].x*RYX + mol[i].y*RYY + mol[i].z*RYZ;
		Q[i + 2 * size] = mol[i].x*RZX + mol[i].y*RZY + mol[i].z*RZZ;
		mol[i].x = Q[i];
		mol[i].y = Q[i + size];
		mol[i].z = Q[i + 2 * size];
	}
}



double MarquesEnantiomers::rmsFitConnect(
	vector<CoordXYZ> &mol1,
	vector<int> &connect1,
	vector<CoordXYZ> &mol2,
	vector<int> &connect2)
{
	double rmsfit = 0.0e0;
	double norm = 0.0e0;
	for (size_t i = 0; i < mol1.size(); i++)
	{
		int i1 = connect1[i];
		int i2 = connect2[i];
		double weigh = 1.0e0; //wfit(i)
		double xr = mol1[i1].x - mol2[i2].x;
		double yr = mol1[i1].y - mol2[i2].y;
		double zr = mol1[i1].z - mol2[i2].z;
		double dist2 = xr*xr + yr*yr + zr*zr;
		norm = norm + weigh;
		double rmsterm = dist2 * weigh;
		rmsfit += rmsterm;
	}
	return sqrt(rmsfit/norm);
}

vector< vector<double> > MarquesEnantiomers::distanceMatrix(
	vector<CoordXYZ> &mol1,
	vector<CoordXYZ> &mol2)
{
	vector< vector<double> > distMatr;
	for (size_t i = 0; i < mol1.size(); i++)
	{
		vector<double> distI(mol1.size());
		for (size_t j = 0; j < mol1.size(); j++)
		{
			double rx = mol1[i].x - mol2[j].x;
			double ry = mol1[i].y - mol2[j].y;
			double rz = mol1[i].z - mol2[j].z;
			distI[j] = sqrt(rx*rx + ry*ry + rz*rz);
			if (mol1[i].atomlabel != mol2[j].atomlabel)
				distI[j] += elementDiff;
		}
		distMatr.push_back(distI);
	}
	return distMatr;
}

void MarquesEnantiomers::quatFit(
	const vector<CoordXYZ> &mol1,
	const vector<int> &connect1,
	vector<CoordXYZ> &mol2,
	const vector<int> &connect2)
{
	int nFit = mol1.size();
	vector< vector<double> > c(4);
	for (size_t i = 0; i < 4; i++)
	{
		c[i].resize(4);
		for (size_t j = 0; j < 4; j++)
			c[i][j] = 0.0e0;
	}

	for (int i = 0; i < nFit; i++)
	{
		int i1 = connect1[i];
		int i2 = connect2[i];
		double weigh = 1.0e0; // possible weight some coordinates
		double xm, ym, zm, xp, yp, zp;
		xm = mol1[i1].x - mol2[i2].x;
		ym = mol1[i1].y - mol2[i2].y;
		zm = mol1[i1].z - mol2[i2].z;
		xp = mol1[i1].x + mol2[i2].x;
		yp = mol1[i1].y + mol2[i2].y;
		zp = mol1[i1].z + mol2[i2].z;

		c[0][0] += weigh*(xm*xm + ym*ym + zm*zm);
		c[0][1] += weigh*(yp*zm - ym*zp);
		c[0][2] += weigh*(xm*zp - xp*zm);
		c[0][3] += weigh*(xp*ym - xm*yp);
		c[1][1] += weigh*(yp*yp + zp*zp + xm*xm);
		c[1][2] += weigh*(xm*ym - xp*yp);
		c[1][3] += weigh*(xm*zm - xp*zp);
		c[2][2] += weigh*(xp*xp + zp*zp + ym*ym);
		c[2][3] += weigh*(ym*zm - yp*zp);
		c[3][3] += weigh*(xp*xp + yp*yp + zm*zm);
	}

	vector<double> d;
	vector< vector<double> > v;
	vector< vector<double> > saveC = c;

	jacobiDiagonalization(c, v, d);
	vector<double> q(4);
	q[0] = v[0][0];
	q[1] = v[1][0];
	q[2] = v[2][0];
	q[3] = v[3][0];

	vector< vector<double> > rot(3);
	rot[0].resize(3);
	rot[1].resize(3);
	rot[2].resize(3);
	
	rot[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	rot[1][0] = 2.0e0 * (q[1] * q[2] - q[0] * q[3]);
	rot[2][0] = 2.0e0 * (q[1] * q[3] + q[0] * q[2]);
	rot[0][1] = 2.0e0 * (q[2] * q[1] + q[0] * q[3]);
	rot[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	rot[2][1] = 2.0e0 * (q[2] * q[3] - q[0] * q[1]);
	rot[0][2] = 2.0e0 * (q[3] * q[1] - q[0] * q[2]);
	rot[1][2] = 2.0e0 * (q[3] * q[2] + q[0] * q[1]);
	rot[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
	//
	//     rotate second molecule to best fit with first(reference) molecule
	//
	for (size_t i = 0; i < mol2.size(); i++)
	{
		double xrot = mol2[i].x * rot[0][0] + mol2[i].y * rot[0][1] + mol2[i].z * rot[0][2];
		double yrot = mol2[i].x * rot[1][0] + mol2[i].y * rot[1][1] + mol2[i].z * rot[1][2];
		double zrot = mol2[i].x * rot[2][0] + mol2[i].y * rot[2][1] + mol2[i].z * rot[2][2];
		mol2[i].x = xrot;
		mol2[i].y = yrot;
		mol2[i].z = zrot;
	}
}

void MarquesEnantiomers::jacobiDiagonalization(
	vector< vector<double> > &a,
	vector< vector<double> > &eigenVectors,
	vector<double> &eigenValues)
{
	int n = a.size();
	int maxrot = 100;
	int nrot = 0;
	vector< vector<double> > v;
	v = a;
	for (int ip = 0; ip < n; ip++)
	{
		for (int iq = 0; iq < n; iq++)
		{
			v[ip][iq] = 0.0e0;
		}
		v[ip][ip] = 1.0e0;
	}
	vector<double> b(n);
	vector<double> d(n);
	vector<double> z(n);
	for (int ip = 0; ip < n; ip++)
	{
		b[ip] = a[ip][ip];
		d[ip] = b[ip];
		z[ip] = 0.0e0;
	}

	double sm, tresh, g, h, t, theta, c, s, tau;
	for (int i = 0; i < maxrot; i++)
	{
		sm = 0.0e0;
		for (int ip = 0; ip < n - 1; ip++)
		{
			for (int iq = ip + 1; iq < n; iq++)
			{
				sm += abs(a[ip][iq]);
			}
		}

		if (sm == 0.0e0)
			break;
		if (i < 3)
			tresh = 0.2e0 * sm / ((double)n * (double)n);
		else
			tresh = 0.0e0;
			
		for (int ip = 0; ip < n - 1; ip++)
		{
			for (int iq = ip + 1; iq < n; iq++)
			{
				g = 100.0e0 * abs(a[ip][iq]);
				if ((i > 3)
					&& ((abs(d[ip]) + g) == abs(d[ip]))
					&& ((abs(d[iq]) + g) == abs(d[iq])))
				{
					a[ip][iq] = 0.0e0;
				}
				else if (abs(a[ip][iq]) > tresh)
				{
					h = d[iq] - d[ip];
					if ((abs(h) + g) == abs(h))
						t = a[ip][iq] / h;
					else
					{
						theta = 0.5e0 * h / a[ip][iq];
						t = 1.0e0 / (abs(theta) + sqrt(1.0e0 + theta * theta));
						if (theta < 0.0e0)
							t = -t;
					}
					c = 1.0e0 / sqrt(1.0e0 + t*t);
					s = t * c;
					tau = s / (1.0e0 + c);
					h = t * a[ip][iq];
					z[ip] = z[ip] - h;
					z[iq] = z[iq] + h;
					d[ip] = d[ip] - h;
					d[iq] = d[iq] + h;
					a[ip][iq] = 0.0e0;
					for (int j = 0; j < ip ; j++)
					{
						g = a[j][ip];
						h = a[j][iq];
						a[j][ip] = g - s*(h + g*tau);
						a[j][iq] = h + s*(g - h*tau);
					}
					for (int j = ip + 1; j < iq ; j++)
					{
						g = a[ip][j];
						h = a[j][iq];
						a[ip][j] = g - s*(h + g*tau);
						a[j][iq] = h + s*(g - h*tau);
					}
					for (int j = iq + 1; j < n; j++)
					{
						g = a[ip][j];
						h = a[iq][j];
						a[ip][j] = g - s*(h + g*tau);
						a[iq][j] = h + s*(g - h*tau);
					}
					for (int j = 0; j < n; j++)
					{
						g = v[j][ip];
						h = v[j][iq];
						v[j][ip] = g - s*(h + g*tau);
						v[j][iq] = h + s*(g - h*tau);
					}
					nrot++;
				}
			}
		}
		for (int ip = 0; ip < n; ip++)
		{
			b[ip] = b[ip] + z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0e0;
		}
	}
	if (nrot == maxrot)
	{
		cout << "JACOBI DIDN'T CONVERGED" << endl;
		exit(1);
	}
	//
	// ORGANIZAR OS AUTOVALORES E AUTOVETORES
	//
	int k;
	double p;
	for (int i = 0; i < n - 1; i++)
	{
		k = i;
		p = d[i];
		for (int j = i + 1; j < n; j++)
		{
			if (d[j] < p)
			{
				k = j;
				p = d[j];
			}
		}
		if (k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for (int j = 0; j < n; j++)
			{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}

	eigenVectors = v;
	eigenValues = d;
}

void MarquesEnantiomers::printCoordXyz(vector<CoordXYZ> &coord)
{
	ofstream xyz_("teste.xyz");
	xyz_ << coord.size() << endl
		<< "titulo" << endl;
	for (size_t i = 0; i < coord.size(); i++)
	{
		xyz_ << "H " << "  "
			<< coord[i].x << "  "
			<< coord[i].y << "  "
			<< coord[i].z << endl;
	}
}


void MarquesEnantiomers::maxMinBetweenTwoNumber(
	const int n,
	const int m,
	int &iMax,
	int &iMin)
{
	if (n == m)
		iMax = iMin = n;
	else if (n > m)
	{
		iMax = n;
		iMin = m;
	}
	else
	{
		iMax = m;
		iMin = n;
	}
}























































////////////////////////////////////////////////////////


double MarquesEnantiomers::setMass(string atomLabel)
{
	if (atomLabel == "H")
		return 1.00794e0;
	else if (atomLabel == "He")
		return 4.002602e0;
	else if (atomLabel == "Li")
		return 6.941e0;
	else if (atomLabel == "Be")
		return 9.01218e0;
	else if (atomLabel == "B")
		return 10.811e0;
	else if (atomLabel == "C")
		return 12.011e0;
	else if (atomLabel == "N")
		return 14.00674e0;
	else if (atomLabel == "O")
		return 15.9994e0;
	else if (atomLabel == "F")
		return 18.998403e0;
	else if (atomLabel == "Ne")
		return 20.1797e0;
	else if (atomLabel == "Na")
		return 22.989768e0;
	else if (atomLabel == "Mg")
		return 24.305e0;
	else if (atomLabel == "Al")
		return 26.981539e0;
	else if (atomLabel == "Si")
		return 28.0855e0;
	else if (atomLabel == "P")
		return 30.973762e0;
	else if (atomLabel == "S")
		return 32.066e0;
	else if (atomLabel == "Cl")
		return 35.4527e0;
	else if (atomLabel == "Ar")
		return 39.948e0;
	else if (atomLabel == "K")
		return 39.0983e0;
	else if (atomLabel == "Ca")
		return 40.078e0;
	else if (atomLabel == "Sc")
		return 44.95591e0;
	else if (atomLabel == "Ti")
		return 47.88e0;
	else if (atomLabel == "V")
		return 50.9415e0;
	else if (atomLabel == "Cr")
		return 51.9961e0;
	else if (atomLabel == "Mn")
		return 54.93805e0;
	else if (atomLabel == "Fe")
		return 55.847e0;
	else if (atomLabel == "Co")
		return 58.9332e0;
	else if (atomLabel == "Ni")
		return 58.6934e0;
	else if (atomLabel == "Cu")
		return 63.546e0;
	else if (atomLabel == "Zn")
		return 65.39e0;
	else if (atomLabel == "Ga")
		return 69.723e0;
	else if (atomLabel == "Ge")
		return 72.61e0;
	else if (atomLabel == "As")
		return 74.92159e0;
	else if (atomLabel == "Se")
		return 78.96e0;
	else if (atomLabel == "Br")
		return 79.904e0;
	else if (atomLabel == "Kr")
		return 83.8e0;
	else if (atomLabel == "Rb")
		return 85.4678e0;
	else if (atomLabel == "Sr")
		return 87.62e0;
	else if (atomLabel == "Y")
		return 88.90585e0;
	else if (atomLabel == "Zr")
		return 91.224e0;
	else if (atomLabel == "Nb")
		return 92.90638e0;
	else if (atomLabel == "Mo")
		return 95.94e0;
	else if (atomLabel == "Tc")
		return 97.9072e0;
	else if (atomLabel == "Ru")
		return 101.07e0;
	else if (atomLabel == "Rh")
		return 102.9055e0;
	else if (atomLabel == "Pd")
		return 106.42e0;
	else if (atomLabel == "Ag")
		return 107.8682e0;
	else if (atomLabel == "Cd")
		return 112.411e0;
	else if (atomLabel == "In")
		return 114.818e0;
	else if (atomLabel == "Sn")
		return 118.71e0;
	else if (atomLabel == "Sb")
		return 121.760e0;
	else if (atomLabel == "Te")
		return 127.6e0;
	else if (atomLabel == "I")
		return 126.90447e0;
	else if (atomLabel == "Xe")
		return 131.29e0;
	else if (atomLabel == "Cs")
		return 132.90543e0;
	else if (atomLabel == "Ba")
		return 137.327e0;
	else if (atomLabel == "La")
		return 138.9055e0;
	else if (atomLabel == "Ce")
		return 140.115e0;
	else if (atomLabel == "Pr")
		return 140.90765e0;
	else if (atomLabel == "Nd")
		return 144.24e0;
	else if (atomLabel == "Pm")
		return 144.9127e0;
	else if (atomLabel == "Sm")
		return 150.36e0;
	else if (atomLabel == "Eu")
		return 151.965e0;
	else if (atomLabel == "Gd")
		return 157.25e0;
	else if (atomLabel == "Tb")
		return 158.92534e0;
	else if (atomLabel == "Dy")
		return 162.50e0;
	else if (atomLabel == "Ho")
		return 164.93032e0;
	else if (atomLabel == "Er")
		return 167.26e0;
	else if (atomLabel == "Tm")
		return 168.93421e0;
	else if (atomLabel == "Yb")
		return 173.04e0;
	else if (atomLabel == "Lu")
		return 174.967e0;
	else if (atomLabel == "Hf")
		return 178.49e0;
	else if (atomLabel == "Ta")
		return 180.9479e0;
	else if (atomLabel == "W")
		return 183.84e0;
	else if (atomLabel == "Re")
		return 186.207e0;
	else if (atomLabel == "Os")
		return 190.23e0;
	else if (atomLabel == "Ir")
		return 192.22e0;
	else if (atomLabel == "Pt")
		return 195.08e0;
	else if (atomLabel == "Au")
		return 196.96654e0;
	else if (atomLabel == "Hg")
		return 200.59e0;
	else if (atomLabel == "Tl")
		return 204.3833e0;
	else if (atomLabel == "Pb")
		return 207.2e0;
	else if (atomLabel == "Bi")
		return 208.98037e0;
	else if (atomLabel == "Po")
		return 208.9824e0;
	else if (atomLabel == "Rn")
		return 209.9871e0;
	else if (atomLabel == "Fr")
		return 222.0176e0;
	else if (atomLabel == "Ra")
		return 223.0197e0;
	else if (atomLabel == "Ac")
		return 226.0254e0;
	else if (atomLabel == "Th")
		return 227.0278e0;
	else if (atomLabel == "Pa")
		return 232.0381e0;
	else if (atomLabel == "U")
		return 231.03588e0;
	else if (atomLabel == "Np")
		return 238.0289e0;
	else if (atomLabel == "Pu")
		return 237.048e0;
	else if (atomLabel == "Am")
		return 244.0642e0;
	else if (atomLabel == "Cm")
		return 243.0614e0;
	else if (atomLabel == "Bk")
		return 247.0703e0;
	else if (atomLabel == "Cf")
		return 249.0703e0;
	else if (atomLabel == "Es")
		return 251.0796e0;
	else if (atomLabel == "Fm")
		return 252.083e0;
	else if (atomLabel == "Md")
		return 257.0951e0;
	else if (atomLabel == "No")
		return 258.1e0;
	else if (atomLabel == "Lr")
		return 259.1009e0;
	else
	{
		cout << "atomLabel:  " << atomLabel << "  not found" << endl;
		exit(1);
	}
}
