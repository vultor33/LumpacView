#include "Combinations.h"

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace std;

Combinations::Combinations(){}

Combinations::~Combinations(){}

void Combinations::doAllCombinations()
{

	vector<string> elem(12);
	elem[0] = "m01";
	elem[1] = "m02";
	elem[2] = "m03";
	elem[3] = "m04";
	elem[4] = "m05";
	elem[5] = "m06";
	elem[6] = "B01";
	elem[7] = "B02";
	elem[8] = "B03";
	elem[9] = "C01";
	elem[10] = "C02";
	elem[11] = "C03";
	vector<int> coord(12);
	coord[0] = 1;
	coord[1] = 1;
	coord[2] = 1;
	coord[3] = 1;
	coord[4] = 1;
	coord[5] = 1;
	coord[6] = 2;
	coord[7] = 2;
	coord[8] = 2;
	coord[9] = 2;
	coord[10] = 2;
	coord[11] = 2;
	int size = 12;
	vector<string> allNames;
	string name = "";
	int actualCoord = 0;
	ofstream arc_("combinacoes-do-6.txt");
	for (int i1 = 0; i1 < size; i1++)
	{
		string nameI1 = elem[i1];
		name = nameI1;
		int coordI1 = coord[i1];
		int totalCoord = coordI1;
		if (totalCoord == 6)
		{
			arc_ << name << endl;
		}
		else if (actualCoord > 6)
		{
			actualCoord++;
		}


		for (int i2 = i1; i2 < size; i2++)
		{
			string nameI2 = elem[i2];
			name = nameI1 + nameI2;
			int coordI2 = coord[i2];
			int totalCoord = coordI1 + coordI2;
			if (totalCoord == 6)
			{
				arc_ << name << endl;
			}
			else if (actualCoord > 6)
			{
				actualCoord++;
			}

			for (int i3 = i2; i3 < size; i3++)
			{
				string nameI3 = elem[i3];
				name = nameI1 + nameI2 + nameI3;
				int coordI3 = coord[i3];
				int totalCoord = coordI1 + coordI2 + coordI3;
				if (totalCoord == 6)
				{
					arc_ << name << endl;
				}
				else if (actualCoord > 6)
				{
					actualCoord++;
				}


				for (int i4 = i3; i4 < size; i4++)
				{
					string nameI4 = elem[i4];
					name = nameI1 + nameI2 + nameI3 + nameI4;
					int coordI4 = coord[i4];
					int totalCoord = coordI1 + coordI2 + coordI3 + coordI4;
					if (totalCoord == 6)
					{
						arc_ << name << endl;
					}
					else if (actualCoord > 6)
					{
						actualCoord++;
					}


					for (int i5 = i4; i5 < size; i5++)
					{
						string nameI5 = elem[i5];
						name = nameI1 + nameI2 + nameI3 + nameI4 + nameI5;
						int coordI5 = coord[i5];
						int totalCoord = coordI1 + coordI2 + coordI3 + coordI4 + coordI5;
						if (totalCoord == 6)
						{
							arc_ << name << endl;
						}
						else if (actualCoord > 6)
						{
							actualCoord++;
						}


						for (int i6 = i5; i6 < size; i6++)
						{
							string nameI6 = elem[i6];
							name = nameI1 + nameI2 + nameI3 + nameI4 + nameI5 + nameI6;
							int coordI6 = coord[i6];
							int totalCoord = coordI1 + coordI2 + coordI3 + coordI4 + coordI5 + coordI6;
							if (totalCoord == 6)
							{
								arc_ << name << endl;
							}
							else if (actualCoord > 6)
							{
								actualCoord++;
							}
						}
					}
				}
			}
		}
	}
	arc_.close();
}


vector< vector<int> > Combinations::stringToNumber(string entryString)
{
	//cortar os codigos
	vector<string> allCodes;
	string code;
	for (size_t i = 0; i < entryString.size(); i++)
	{
		if (((i % 3) == 0) && (i != 0))
		{
			allCodes.push_back(code);
			cout << code << endl;
			code = "";
		}
		code += entryString[i];
	}
	allCodes.push_back(code);
	cout << code << endl;

	vector<int> typeCode1;
	vector<int> typeCode2;
	vector<int> typeCode3;
	addDifferent(
		codeToType(allCodes[0]),
		typeCode1,
		typeCode2,
		typeCode3);
	for (size_t i = 1; i < allCodes.size(); i++)
	{
		if (allCodes[i] == allCodes[i - 1])
		{
			addEqual(
				codeToType(allCodes[i]),
				typeCode1,
				typeCode2,
				typeCode3);
		}
		else
		{
			addDifferent(
				codeToType(allCodes[i]),
				typeCode1,
				typeCode2,
				typeCode3);
		}
	}

	sort(typeCode1.begin(), typeCode1.end());
	sort(typeCode2.begin(), typeCode2.end());
	sort(typeCode3.begin(), typeCode3.end());

	vector< vector<int> > types;
	types.push_back(typeCode1);
	types.push_back(typeCode2);
	types.push_back(typeCode3);

	// parece que esta ok
	return types;
}

int Combinations::codeToType(string code)
{
	if (code[0] == 'm')
		return 0;
	if (code[0] == 'B')
		return 1;
	if (code[0] == 'C')
		return 2;
}

void Combinations::addEqual(
	int codeNumber,
	std::vector<int> & typeCode1,
	std::vector<int> & typeCode2,
	std::vector<int> & typeCode3
	)
{
	if (codeNumber == 0)
	{
		typeCode1[typeCode1.size() - 1]++;
	}
	else if (codeNumber == 1)
	{
		typeCode2[typeCode2.size() - 1]++;
	}
	else if (codeNumber == 2)
	{
		typeCode3[typeCode3.size() - 1]++;
	}
}

void Combinations::addDifferent(
	int codeNumber,
	std::vector<int> & typeCode1,
	std::vector<int> & typeCode2,
	std::vector<int> & typeCode3
	)
{
	int newCode = 1;
	if (codeNumber == 0)
	{
		typeCode1.push_back(newCode);
	}
	else if (codeNumber == 1)
	{
		typeCode2.push_back(newCode);
	}
	else if (codeNumber == 2)
	{
		typeCode3.push_back(newCode);
	}
}

