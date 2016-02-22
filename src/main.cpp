#include <stdlib.h>
#include <iostream>

#include "MyExceptions.h"
#include "Coordstructs.h"
#include "ReadInput.h"

using namespace std;

int main()
{
	ReadInput readInp_;
	try	{
		readInp_.readLumpacViewInput();
	} catch (MyExceptions& e)	{
		cout << e.what() << endl;
		return 1;
	}
	catch (...) {
		cout << "unknown error on input - check it or contact developers" << endl;
		return 1;
	}

#ifdef _DEBUG
	readInp_.rePrintInput();
#endif
	// passar as informacoes atraves de referencias.


	// aplicar a matematica do simas.
	// no caso do monodentdo o centro de massa funciona bem
	bool sucess;
	bool terminateExecution = false;
	for (size_t i = 0; i < readInp_.allLigands.size(); i++)
	{
		sucess = readInp_.allLigands[i].initializeLigand();
		if (!sucess) 
		{
			cout << "Can't build ligand:  " << i << " check input" << endl;
			terminateExecution = true;
		}
	}
	if (terminateExecution)
		return 1;





	cin.get();
	return 0;
}



// Estou usando o preprocessor _DEBUG pra testes,
// só vendo se a parada funciona, qualquer e so dar um find e apagar

/*
ERROR HANDLING
- Quando o destructor nao e chamado:
- Destructor of the class is not called if exception is thrown in its constructor.
- Exception is automatically re-thrown if caught in construction initialization list catch block.
*/

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



