#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
	stringstream ssReadInput;
	ssReadInput << argv[1] << "  " << argv[2] << "  " << argv[3];
	string responseName;
	int kInit, kFinal;
	ssReadInput >> responseName >> kInit >> kFinal;

	ifstream resp_(responseName.c_str());
	string line;
	int k = 1;
	ofstream roda_("roda.x");
	roda_ << "#!/bin/bash" << endl;
	while(getline(resp_,line))
	{
		if((k >= kInit) && (k <= kFinal))
		{
			stringstream convert;
			convert << line;
			string name;
			convert >> name;
			roda_ << "./script-run-flag-pc " << name << endl;

		} 
		k++;
	}
	system("chmod u+x roda.x");
	roda_.close();
	resp_.close();

	return 0;

}




