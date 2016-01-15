#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <time.h>

#include "Coordstructs.h"
#include "ReadInput.h"

using namespace std;

int main()
{
	ReadInput readInp_;

	readInp_.rePrintInput();

	cin.get();
	return 0;
}


/*
to string
std::string to_string(const int& n);
string ReadInput::to_string( const int& n )
{
std::ostringstream stm ;
stm << n ;
return stm.str() ;
}



for (int i = 0; i < (int)ri_.numberLigands.size(); i++)
{
for (int j = 0; j < ri_.numberLigands[i]; j++)
{
Ligand auxOptLigand; //com 8 coordXYZ
for (int l = 0; l < (int)initialConfiguration[m].coord.size(); l++) //ligand i -> 8 atoms
{
m++;
CoordXYZ auxCoord;
auxCoord = optimizedMol_.coord[k];
k++;
auxCoord.x -= centerX;
auxCoord.y -= centerY;
auxCoord.z -= centerZ;
auxOptLigand.coord.push_back(auxCoord);
}
allOptimized.push_back(auxOptLigand); // add ligand i
}
}











for (int i = 0; i < 50; i++)
{

	Ligand h2o = ri_.readConfigurations("h2o-ligante");
	Ligand h2o2 = h2o;
	Ligand h2o3 = h2o;
	Ligand h2o4 = h2o;
	Ligand h2o5 = h2o;
	Ligand h2o6 = h2o;
	Ligand h2o7 = h2o;
	Ligand h2o8 = h2o;

	PositionChange pos_;
	pos_.initialPosition(h2o);
	pos_.initialPosition(h2o2);
	pos_.initialPosition(h2o3);
	pos_.initialPosition(h2o4);
	pos_.initialPosition(h2o5);
	pos_.initialPosition(h2o6);
	pos_.initialPosition(h2o7);
	pos_.initialPosition(h2o8);

	vector<Ligand> all(8);
	all[0] = h2o;
	all[1] = h2o2;
	all[2] = h2o3;
	all[3] = h2o4;
	all[4] = h2o5;
	all[5] = h2o6;
	all[6] = h2o7;
	all[7] = h2o8;

	pr_.printAllMol(all);
}
*/



