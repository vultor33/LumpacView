#include "MyExceptions.h"

#include <string>

using namespace std;

MyExceptions::MyExceptions(int iError_in)
{
	iError = iError_in;
}

MyExceptions::MyExceptions(string customMessage_in)
{
	iError = 0;
	customMessage = customMessage_in;
}


MyExceptions::~MyExceptions(){}

const char * MyExceptions::what() const throw()
{
	{
		switch (iError)
		{
		case 0:
			return customMessage.c_str();
		case 1:
			return "Error 1 - Ligand file not found";
		case 2:
			return "Error 2 - 000"; //avaible
		case 3:
			return "Error 3 - Wrong number of atoms at ligand file";
		case 4:
			return "Error 4 - LumpacViewInput.txt not found";

		default:
			return "unknow error - contact developers";
		}
	}
}

