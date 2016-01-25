#include "MyExceptions.h"

using namespace std;

MyExceptions::MyExceptions(int iError_in)
	:iError(iError_in)
{
}

MyExceptions::~MyExceptions(){}

const char * MyExceptions::what() const throw()
{
	{
		switch (iError)
		{
		case 1:
			return "Error 1 - Ligand file not found";
		case 2:
			return "Error 2 - Problem on ligands, check input";
		case 3:
			return "Error 3 - Wrong number of atoms at ligand file";

		default:
			return "unknow error - contact developers";
		}
	}
}

