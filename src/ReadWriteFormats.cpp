#include "ReadWriteFormats.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

ReadWriteFormats::~ReadWriteFormats(){}

// number / bidentares   --> format
void ReadWriteFormats::readAtomTypesAndBidentateChosenFile(
	ifstream & file_,
	vector<int> & atomTypes,
	vector<int> & bidentateChosen,
	int systemSize,
	int nBidentates)
{
	string line;
	getline(file_, line);
	stringstream auxline;
	auxline << line;
	for (int i = 0; i < systemSize; i++)
	{
		int aux;
		auxline >> aux;
		atomTypes.push_back(aux - 1);
	}
	if (nBidentates != 0)
	{
		string dummy;
		auxline >> dummy;
		for (int i = 0; i < 2 * nBidentates; i++)
		{
			int aux;
			auxline >> aux;
			bidentateChosen.push_back(aux - 1);
		}
	}
}


// {SAPR-8 [Ma2b2c2] [1 2 3 4 5 6 7 8] Aa} --> format
vector<int> ReadWriteFormats::readCauchyNotationsEnantiomers(
	ifstream & openendFile_,
	int size)
{
	vector<int> notation;
	if (openendFile_.eof())
		return notation;

	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;
	notation.resize(size);

	size_t brack1Temp = auxline.find("]");
	size_t brack1 = auxline.find("[", brack1Temp + 1, 1);
	size_t brack2 = auxline.find("]", brack1Temp + 1, 1);
	string permString = auxline.substr(brack1 + 1, brack2 - brack1 - 1);
	stringstream line;
	line << permString;
	for (size_t i = 0; i < size; i++)
	{
		line >> notation[i];
		//notation[i]--; --- observar isso aqui
	}
	return notation;
}

// 3 vec (mono, bsym and bass) --> m01m02m03m04m06B01C01 
string ReadWriteFormats::codeToString(std::vector < std::vector<int> > & codeLine)
{
	string name = "";
	for (size_t i = 0; i < codeLine.size(); i++)
	{
		int startCount;
		if (i == 0)
			startCount = 0;
		else if (i == 1)
			startCount = 12;
		else if (i == 2)
			startCount = 18;
		for (size_t j = 0; j < codeLine[i].size(); j++)
		{
			for (int k = 0; k < codeLine[i][j]; k++)
			{
				name += elem[j + startCount];
			}
		}
	}
	return name;
}

// 3 vec (mono, bsym and bass) --> Ma2bcd(AA)(AB)
string ReadWriteFormats::newCodeToString(vector < vector<int> > & codeLine)
{
	string name = "";
	for (size_t i = 0; i < codeLine.size(); i++)
	{
		int k = 0;
		int startCount;
		if (i == 0)
			startCount = 0;
		else if (i == 1)
			startCount = 12;
		else if (i == 2)
			startCount = 18;
		for (int j = (int)codeLine[i].size() - 1; j > -1; j--)
		{
			stringstream convert;
			if (codeLine[i][j] > 1)
			{
				convert << elemNew[k + startCount]
					<< codeLine[i][j];
				k++;
			}
			else
			{
				convert << elemNew[k + startCount];
				k++;
			}
			name += convert.str();
		}
	}
	return name;
}


// (m, B e C)  versao do allMolecular
vector< vector<int> > ReadWriteFormats::stringToNumber(string entryString)
{
	vector<string> allCodes;
	string code;
	for (size_t i = 0; i < entryString.size(); i++)
	{
		if (((i % 3) == 0) && (i != 0))
		{
			allCodes.push_back(code);
			code = "";
		}
		code += entryString[i];
	}
	allCodes.push_back(code);

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

	return types;
}

int ReadWriteFormats::codeToType(string code)
{
	if (code[0] == 'm')
		return 0;
	else if (code[0] == 'B')
		return 1;
	else if (code[0] == 'C')
		return 2;

	return -1;
}

void ReadWriteFormats::addEqual(
	int codeNumber,
	std::vector<int> & typeCode1,
	std::vector<int> & typeCode2,
	std::vector<int> & typeCode3)
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

void ReadWriteFormats::addDifferent(
	int codeNumber,
	std::vector<int> & typeCode1,
	std::vector<int> & typeCode2,
	std::vector<int> & typeCode3)
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

// {a (AA) (AB)}  versao do isomers to mol
int ReadWriteFormats::CompositionToNumbers(
	string entryString, 
	int & nBidentates)
{
	// cut ligands codes
	vector<string> allCodes;
	int coordination = 0;
	nBidentates = 0;
	bool bidentate;
	for (size_t i = 1; i < entryString.size(); i++)
	{
		char codeI = entryString[i];
		bidentate = false;
		if (codeI == '(')
		{
			i += 4;
			bidentate = true;
		}
		else
			i += 1;
		if ((i == entryString.size()) ||
			(!isdigit(entryString[i])))
		{
			i--;
			if (bidentate)
			{
				nBidentates += 1;
				coordination += 2;
			}
			else
				coordination++;
		}
		else
		{
			stringstream convertI;
			convertI << entryString[i];
			int nTimes;
			convertI >> nTimes;
			if (bidentate)
			{
				nBidentates += nTimes;
				coordination += 2 * nTimes;
			}
			else
			{
				coordination += nTimes;
			}
		}
	}
	return coordination;
}






ReadWriteFormats::ReadWriteFormats()
{
	elem.resize(24);
	elem[0] = "m01";
	elem[1] = "m02";
	elem[2] = "m03";
	elem[3] = "m04";
	elem[4] = "m05";
	elem[5] = "m06";
	elem[6] = "m07";
	elem[7] = "m08";
	elem[8] = "m09";
	elem[9] = "m10";
	elem[10] = "m11";
	elem[11] = "m12";
	elem[12] = "B01";
	elem[13] = "B02";
	elem[14] = "B03";
	elem[15] = "B04";
	elem[16] = "B05";
	elem[17] = "B06";
	elem[18] = "C01";
	elem[19] = "C02";
	elem[20] = "C03";
	elem[21] = "C04";
	elem[22] = "C05";
	elem[23] = "C06";
	elemNew.resize(24);
	elemNew[0] = "a";
	elemNew[1] = "b";
	elemNew[2] = "c";
	elemNew[3] = "d";
	elemNew[4] = "e";
	elemNew[5] = "f";
	elemNew[6] = "g";
	elemNew[7] = "h";
	elemNew[8] = "i";
	elemNew[9] = "j";
	elemNew[10] = "k";
	elemNew[11] = "l";
	elemNew[12] = "(AA)";
	elemNew[13] = "(BB)";
	elemNew[14] = "(CC)";
	elemNew[15] = "(DD)";
	elemNew[16] = "(EE)";
	elemNew[17] = "(FF)";
	elemNew[18] = "(AB)";
	elemNew[19] = "(CD)";
	elemNew[20] = "(EF)";
	elemNew[21] = "(GH)";
	elemNew[22] = "(IJ)";
	elemNew[23] = "(KL)";
}











