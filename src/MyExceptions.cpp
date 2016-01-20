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

		default:
			return "unknow error - contact developers";
		}
	}
}

