#include "UtilityRun.h"

#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

UtilityRun::UtilityRun(){}

UtilityRun::~UtilityRun(){}

void UtilityRun::renameAtomTypes(string responseName)
{
	ifstream response_(responseName.c_str());
	string line;
	while (!response_.eof())
	{
		getline(response_, line);
		if (line == "")
			break;

		stringstream convert;
		convert << line;
		string comp;
		convert >> comp;
		string oldName = comp + "---atomTypes.txt";
		string newName = "final-" + comp + "-atomTypes";
		int result = rename(oldName.c_str(), newName.c_str());
		if (result == 0)
		{
			cout << comp << " ...  ...  ...  Ok" << endl;
		}
		else
		{
			cout << "Erro na hora de renomear" << endl;
			exit(1);
		}
	}
}



