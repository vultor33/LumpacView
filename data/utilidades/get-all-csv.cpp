#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

int main()
{

	string skeletonName = "isomers-6";
	int size = 6;
	system("cat *.csv > allIsomers.csv");
	string csvName = "allIsomers";

	ifstream csvFile1_((csvName + ".csv").c_str());
	ifstream skeleton_(skeletonName.c_str());
	ofstream csvFileNumber_((csvName + "-number.csv").c_str());

	string line1;
	int k1 = 1;
	while(getline(csvFile1_,line1))
	{
		if(line1 == "")
		{
			csvFileNumber_ << endl;
			continue;
		}


		string aux1, aux2;
		vector<int> permutCsv(size);
		stringstream lineCsv;
		lineCsv << line1;
		lineCsv >> aux1 >> aux2;
		for(int i = 0; i < size; i++)
			lineCsv >> permutCsv[i];

		string line2;
		while(getline(skeleton_,line2))
		{
			if(line2 == "")
				continue;

			vector<int> permutSkeleton(size);
			stringstream lineSkel;
			lineSkel << line2;
			for(int i = 0; i < size; i++)
				lineSkel >> permutSkeleton[i];
			
			if(permutSkeleton == permutCsv)
			{
				csvFileNumber_ << aux1 << "  " 
				<< aux2 << " " << k1 << " ; ";
				for(int i = 0; i < size; i++)
					csvFileNumber_ << permutCsv[i] << "  ";

				csvFileNumber_ << endl;
				k1++;
				break;
			}
			k1++;
		}
	}
	csvFileNumber_.close();
	skeleton_.close();
	csvFile1_.close();

	system("rm allWeights.txt");
	system("cat *weights* > allWeights.txt");
	//system("rm *weights*");
	ifstream weightsFile_("allWeights.txt");
	string line;
	vector<int> allWeights;
	while(getline(weightsFile_,line))
		{
		stringstream readLine;
		readLine << line;
		int number;
		readLine >> number;
		allWeights.push_back(number);
	}
	sort(allWeights.begin(),allWeights.end());


	ifstream csvNumber_((csvName + "-number.csv").c_str());
	ofstream csvFinal_((csvName + "-final.csv").c_str());
	int kAll = 0;
	while(getline(csvNumber_,line))
	{
		if(line == "")
		{
			csvFinal_ << endl;
			continue;
		}

		string dum1, dum2;
		vector<int> permutCsv(size);
		int currentWeight, permutIndex;
		stringstream csvLine;
		csvLine << line;
		csvLine >> currentWeight >> dum1 >> permutIndex >> dum2;
		for(int i = 0; i < size; i++)
			csvLine >> permutCsv[i];

		while(permutIndex == allWeights[kAll])
		{
			currentWeight++;
			kAll++;
		}

		csvFinal_ << currentWeight << "  ;  " << permutIndex << "  ;  ";
		for(int i = 0; i < size; i++)
			csvFinal_ << permutCsv[i] << "  ";
		csvFinal_ << endl;

	}

	return 0;

}


