#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <unistd.h>
#include <string>
#include <sstream>

using namespace std;

int main()
{
	while(true)
	{
       	 	sleep(0.1);
        	system("squeue > quee.txt");
        	ifstream in_("quee.txt");
        	string line;
		for(int i = 0; i < 110; i++)
		{
        		getline(in_,line);
		}
		in_.close();
		remove("quee.txt");

        	if(line == "")
			break;
	}
	return 0;
}



