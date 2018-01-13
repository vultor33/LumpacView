#include "AllMolecularFormulas.h"

#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <sstream>

#include "FindIsomers.h"
#include "ReadWriteFormats.h"

using namespace std;

AllMolecularFormulas::AllMolecularFormulas()
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
	coord.resize(24);
	coord[0] = 1;
	coord[1] = 1;
	coord[2] = 1;
	coord[3] = 1;
	coord[4] = 1;
	coord[5] = 1;
	coord[6] = 1;
	coord[7] = 1;
	coord[8] = 1;
	coord[9] = 1;
	coord[10] = 1;
	coord[11] = 1;
	coord[12] = 2;
	coord[13] = 2;
	coord[14] = 2;
	coord[15] = 2;
	coord[16] = 2;
	coord[17] = 2;
	coord[18] = 2;
	coord[19] = 2;
	coord[20] = 2;
	coord[21] = 2;
	coord[22] = 2;
	coord[23] = 2;

	allLigandsNames.resize(36);
	allLigandsNames[0] = "auxLigands/Lumpac-View-Dummy-Ligand-M1";
	allLigandsNames[1] = "auxLigands/Lumpac-View-Dummy-Ligand-M2";
	allLigandsNames[2] = "auxLigands/Lumpac-View-Dummy-Ligand-M3";
	allLigandsNames[3] = "auxLigands/Lumpac-View-Dummy-Ligand-M4";
	allLigandsNames[4] = "auxLigands/Lumpac-View-Dummy-Ligand-M5";
	allLigandsNames[5] = "auxLigands/Lumpac-View-Dummy-Ligand-M6";
	allLigandsNames[6] = "auxLigands/Lumpac-View-Dummy-Ligand-M7";
	allLigandsNames[7] = "auxLigands/Lumpac-View-Dummy-Ligand-M8";
	allLigandsNames[8] = "auxLigands/Lumpac-View-Dummy-Ligand-M9";
	allLigandsNames[9] = "auxLigands/Lumpac-View-Dummy-Ligand-M10";
	allLigandsNames[10] = "auxLigands/Lumpac-View-Dummy-Ligand-M11";
	allLigandsNames[11] = "auxLigands/Lumpac-View-Dummy-Ligand-M12";
	allLigandsNames[12] = "auxLigands/Lumpac-View-Dummy-Ligand-B11";
	allLigandsNames[13] = "auxLigands/Lumpac-View-Dummy-Ligand-B21";
	allLigandsNames[14] = "auxLigands/Lumpac-View-Dummy-Ligand-B31";
	allLigandsNames[15] = "auxLigands/Lumpac-View-Dummy-Ligand-B41";
	allLigandsNames[16] = "auxLigands/Lumpac-View-Dummy-Ligand-B51";
	allLigandsNames[17] = "auxLigands/Lumpac-View-Dummy-Ligand-B61";
	allLigandsNames[18] = "auxLigands/Lumpac-View-Dummy-Ligand-B12";
	allLigandsNames[19] = "auxLigands/Lumpac-View-Dummy-Ligand-B22";
	allLigandsNames[20] = "auxLigands/Lumpac-View-Dummy-Ligand-B32";
	allLigandsNames[21] = "auxLigands/Lumpac-View-Dummy-Ligand-B42";
	allLigandsNames[22] = "auxLigands/Lumpac-View-Dummy-Ligand-B52";
	allLigandsNames[23] = "auxLigands/Lumpac-View-Dummy-Ligand-B62";
	allLigandsNames[24] = "auxLigands/Lumpac-View-Dummy-Ligand-C11";
	allLigandsNames[25] = "auxLigands/Lumpac-View-Dummy-Ligand-C21";
	allLigandsNames[26] = "auxLigands/Lumpac-View-Dummy-Ligand-C31";
	allLigandsNames[27] = "auxLigands/Lumpac-View-Dummy-Ligand-C41";
	allLigandsNames[28] = "auxLigands/Lumpac-View-Dummy-Ligand-C51";
	allLigandsNames[29] = "auxLigands/Lumpac-View-Dummy-Ligand-C61";
	allLigandsNames[30] = "auxLigands/Lumpac-View-Dummy-Ligand-C12";
	allLigandsNames[31] = "auxLigands/Lumpac-View-Dummy-Ligand-C22";
	allLigandsNames[32] = "auxLigands/Lumpac-View-Dummy-Ligand-C32";
	allLigandsNames[33] = "auxLigands/Lumpac-View-Dummy-Ligand-C42";
	allLigandsNames[34] = "auxLigands/Lumpac-View-Dummy-Ligand-C52";
	allLigandsNames[35] = "auxLigands/Lumpac-View-Dummy-Ligand-C62";
	ligandsSeparationSize = 6;
}

AllMolecularFormulas::~AllMolecularFormulas(){}

void AllMolecularFormulas::findAllIsomersOnCombinations(string combinationFile)
{
	ofstream printAllCsv_("CombinationsCsv.csv");
	ifstream combFile_(combinationFile.c_str());
	string auxline;
	while (getline(combFile_, auxline))
	{
		if ((auxline == "") || (auxline == "end"))
			break;

		stringstream line;
		line << auxline;
		string code;
		line >> code;
		printInputWithCode(code);
		FindIsomers findIso_;
		int counter = findIso_.start();
		printAllCsv_ << code << " ; " << counter << endl;
	}
	combFile_.close();
	printAllCsv_.close();
}

void AllMolecularFormulas::printInputWithCode(std::string code)
{
	vector<string> ligandNames;
	vector<int> denticity;
	codeToLigands(code, ligandNames, denticity);
	string inputName = "LumpacViewInput.txt";
	remove(inputName.c_str());
	ofstream pInput_(inputName.c_str());
	pInput_ << code + ".xyz" << endl;
	pInput_ << "Eu" << endl << "Eu_spk" << endl;
	pInput_ << ligandNames.size() << endl;
	for (size_t i = 0; i < ligandNames.size(); i++)
		pInput_ << ligandNames[i] << endl;
	pInput_ << angleBidentateCut(ligandNames.size()) << endl;
	pInput_ << denticity.size() << endl;
	for (size_t i = 0; i < denticity.size(); i++)
		pInput_ << denticity[i] << endl;
	pInput_.close();

/*		
teste.xyz
Eu
Eu_spk
5
auxLigands/Lumpac-View-Dummy-Ligand-M1
auxLigands/Lumpac-View-Dummy-Ligand-M1
auxLigands/Lumpac-View-Dummy-Ligand-M2
auxLigands/Lumpac-View-Dummy-Ligand-B11
auxLigands/Lumpac-View-Dummy-Ligand-B12
120
2
3
4
*/
}



void AllMolecularFormulas::doAllCombinations(int nCoordination)
{
	printCoord = nCoordination;
	stringstream convert;
	convert << printCoord;
	string coordName;
	convert >> coordName;
	string fileName = "combinations" + coordName + ".txt";
	maxSize = coord.size();
	vector<string> name(maxSize / 2);
	vector<int> totalCoord(maxSize / 2);
	printCombinations_.open(fileName.c_str());
	for (int i1 = 0; i1 < maxSize; i1++)
	{
		sumUpNameAndPrint(name, totalCoord, 0, elem[i1], coord[i1]);
		for (int i2 = i1; i2 < maxSize; i2++)
		{
			if (printCoord < 2) break;
			sumUpNameAndPrint(name, totalCoord, 1, elem[i2], coord[i2]);
			for (int i3 = i2; i3 < maxSize; i3++)
			{
				if (printCoord < 3) break;
				sumUpNameAndPrint(name, totalCoord, 2, elem[i3], coord[i3]);
				for (int i4 = i3; i4 < maxSize; i4++)
				{
					if (printCoord < 4) break;
					sumUpNameAndPrint(name, totalCoord, 3, elem[i4], coord[i4]);
					for (int i5 = i4; i5 < maxSize; i5++)
					{
						if (printCoord < 5) break;
						sumUpNameAndPrint(name, totalCoord, 4, elem[i5], coord[i5]);
						for (int i6 = i5; i6 < maxSize; i6++)
						{
							if (printCoord < 6) break;
							sumUpNameAndPrint(name, totalCoord, 5, elem[i6], coord[i6]);
							for (int i7 = i6; i7 < maxSize; i7++)
							{
								if (printCoord < 7) break;
								sumUpNameAndPrint(name, totalCoord, 6, elem[i7], coord[i7]);
								for (int i8 = i7; i8 < maxSize; i8++)
								{
									if (printCoord < 8) break;
									sumUpNameAndPrint(name, totalCoord, 7, elem[i8], coord[i8]);
									for (int i9 = i8; i9 < maxSize; i9++)
									{
										if (printCoord < 9) break;
										sumUpNameAndPrint(name, totalCoord, 8, elem[i9], coord[i9]);
										for (int i10 = i9; i10 < maxSize; i10++)
										{
											if (printCoord < 10) break;
											sumUpNameAndPrint(name, totalCoord, 9, elem[i10], coord[i10]);
											for (int i11 = i10; i11 < maxSize; i11++)
											{
												if (printCoord < 11) break;
												sumUpNameAndPrint(name, totalCoord, 10, elem[i11], coord[i11]);
												for (int i12 = i11; i12 < maxSize; i12++)
												{
													if (printCoord < 12) break;
													sumUpNameAndPrint(name, totalCoord, 11, elem[i12], coord[i12]);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	printCombinations_.close();

	clearEqualCombinations(fileName);
	
	remove(fileName.c_str());
}

void AllMolecularFormulas::codeToLigands(
	std::string code,
	std::vector< std::string > & ligandNames,
	std::vector<int> & denticity)
{
	ReadWriteFormats rwf_;
	vector< vector<int> > codeLine = rwf_.compositionToNumberOld(code);
	int l = 0;
	int m = 0;
	for (size_t i = 0; i < codeLine.size(); i++)
	{
		int startCount;
		if (i == 0)
			startCount = 0;
		else if (i == 1)
			startCount = 12;
		else if (i == 2)
			startCount = 24;

		for (size_t j = 0; j < codeLine[i].size(); j++)
		{
			for (int k = 0; k < codeLine[i][j]; k++)
			{
				if (startCount == 0)
				{
					ligandNames.push_back(allLigandsNames[j]);
					l++;
					m++;
				}
				else if (startCount > 0)
				{
					ligandNames.push_back(allLigandsNames[j + startCount]);
					ligandNames.push_back(allLigandsNames[j + ligandsSeparationSize + startCount]);
					denticity.push_back(l);
					denticity.push_back(l + 1);
					l += 2;
					m++;
				}
			}
		}
	}
}

void AllMolecularFormulas::clearEqualCombinations(string combFile)
{
	ifstream combFile_(combFile.c_str());

	vector< vector< vector<int> > > allCodes;

	string auxline;
	getline(combFile_, auxline);
	ReadWriteFormats rwf_;
	vector< vector<int> > code0 = rwf_.compositionToNumberOld(auxline);
	allCodes.push_back(code0);
	while (getline(combFile_, auxline))
	{
		if (auxline == "") break;

		vector< vector<int> > codeLine = rwf_.compositionToNumberOld(auxline);

		if (compareToAll(allCodes, codeLine))
			allCodes.push_back(codeLine);
	}
	combFile_.close();

	string outFileName = "response-" + combFile;
	ofstream out_(outFileName.c_str());

	for (size_t i = 0; i < allCodes.size(); i++)
		out_ << rwf_.codeToString(allCodes[i]) << endl;

	out_.close();
}

/*
string AllMolecularFormulas::codeToString(vector < vector<int> > & codeLine)
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

string AllMolecularFormulas::newCodeToString(vector < vector<int> > & codeLine)
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
*/

bool AllMolecularFormulas::compareToAll(std::vector< std::vector< std::vector<int> > > & allCodes, std::vector< std::vector<int> > &actualCodes)
{
	for (size_t i = 0; i < allCodes.size(); i++)//run codes
	{
		bool compare = true;
		for (size_t j = 0; j < allCodes[i].size(); j++)
		{
			if (actualCodes[j].size() != allCodes[i][j].size())
			{
				compare = false;
				break;
			}
			for (size_t k = 0; k < allCodes[i][j].size(); k++)
			{
				if (actualCodes[j][k] != allCodes[i][j][k])
				{
					compare = false;
					break;
				}
			}
			if (!compare) break;
		}
		if (compare)
			return false;
	}
	return true;
}

void AllMolecularFormulas::sumUpNameAndPrint(
	std::vector< std::string > & name,
	std::vector< int > & totalCoord,
	int position,
	std::string newLigand,
	int coordination)
{
	name[position] = newLigand;
	totalCoord[position] = coordination;
	for (size_t i = position + 1; i < name.size(); i++)
	{
		name[i] = "";
		totalCoord[i] = 0;
	}
	int sumCoord = 0;
	for (size_t i = 0; i < totalCoord.size(); i++)
		sumCoord += totalCoord[i];

	if (sumCoord == printCoord)
	{
		string totalName = "";
		for (size_t i = 0; i < name.size(); i++)
			totalName += name[i];
		printCombinations_ << totalName << endl;
	}
}


/*fredapagar
vector< vector<int> > AllMolecularFormulas::stringToNumber(string entryString)
{
	//cortar os codigos
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

	// parece que esta ok
	// obter o vector de vector para cada um e guardar
	// ler todos os strings e comparar os novos com os anteriores
	// guardar so os diferentes.
	return types;
}

int AllMolecularFormulas::codeToType(string code)
{
	if (code[0] == 'm')
		return 0;
	else if (code[0] == 'B')
		return 1;
	else if (code[0] == 'C')
		return 2;

	return -1;
}

void AllMolecularFormulas::addEqual(
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

void AllMolecularFormulas::addDifferent(
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
*/


int AllMolecularFormulas::angleBidentateCut(size_t coordination)
{
	cout << "On AllMolecularFormulas::angleBidentateCut need to have caution" << endl;
	exit(1);

	switch (coordination)
	{
	case 4:
		return 120;
	case 5:
		return 130;
	case 6:
		return 100;
	case 7:
		return 90;
	case 8:
		return 90;
	case 9:
		return 90;
	case 10:
		return 90;
	case 11:
		return 90;
	case 12:
		return 90;
	default:
		break;
	}
	cout << "Error on AllMolecularFormulas::angleBidentateCut" << endl;
	cout << "Coordination not found" << endl;
	exit(1);

}

