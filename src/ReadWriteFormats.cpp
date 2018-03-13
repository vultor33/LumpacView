#include "ReadWriteFormats.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

#include "AuxMath.h"

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

// labels / bidentares   --> format
void ReadWriteFormats::readAtomTypesAndBidentateChosenFileWithLabels(
	ifstream & file_,
	vector<int> & atomTypes,
	vector<int> & bidentateChosen,
	int systemSize,
	int nBidentates)
{
	string line;
	getline(file_, line);
	stringstream auxline;
	vector<string> atomLabels;
	auxline << line;
	for (int i = 0; i < systemSize; i++)
	{
		string aux;
		auxline >> aux;
		atomLabels.push_back(aux);
	}

	atomTypes = atomStringToAtomTypes(atomLabels);

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
	int rcw;
	string vGroup, code;
	return readCauchyNotation(openendFile_, size, vGroup, rcw, code);


	/*
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
	*/
}

std::vector<int> ReadWriteFormats::readCauchyNotation(
	std::ifstream & openendFile_,
	int size,
	string & vultorGroup,
	int & rcw,
	string & allCode)
{
	vector<int> notation;

	if (openendFile_.eof())
		return notation;

	string auxline;
	getline(openendFile_, auxline);
	if (auxline == "")
		return notation;

	notation.resize(size);

	stringstream convert;
	convert << auxline;
	string dummy;
	convert >> rcw >> dummy;

	size_t brack1Temp = auxline.find("]");
	size_t brack1 = auxline.find("[", brack1Temp + 1, 1);
	size_t brack2 = auxline.find("]", brack1Temp + 1, 1);
	string permString = auxline.substr(brack1 + 1, brack2 - brack1 - 1);
	vultorGroup = auxline.substr(brack2 + 2, 2);
	size_t key1 = auxline.find("{");
	allCode = auxline.substr(key1, allCode.size() - key1);

	stringstream line;
	line << permString;
	for (size_t i = 0; i < size; i++)
	{
		line >> notation[i];
		notation[i]--;
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
vector< vector<int> > ReadWriteFormats::compositionToNumberOld(string entryString)
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
int ReadWriteFormats::compositionToNumbers(
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


string ReadWriteFormats::typeLineToLetters(
	string typeLines,
	std::vector<int> & atomTypes,
	std::vector<int> & bidentateChosen)
{
	vector<int> numberCodes;
	stringstream bidLine;
	int bidNumber = 0;
	for (size_t i = 0; i < typeLines.size(); i++)
	{
		if (typeLines[i] == '/')
		{
			for (size_t j = i + 1; j < typeLines.size(); j++)
			{
				if (typeLines[j] != ' ')
				{
					bidLine << typeLines[j] << " ";
					bidNumber++;
				}
			}
			break;
		}
		if (typeLines[i] != ' ')
		{
			stringstream convert;
			convert << typeLines[i];
			int auxNumberCodes;
			convert >> auxNumberCodes;
			auxNumberCodes--;
			numberCodes.push_back(auxNumberCodes);
		}
	}
	for (int i = 0; i < bidNumber; i++)
	{
		int auxBid;
		bidLine >> auxBid;
		auxBid--;
		bidentateChosen.push_back(auxBid);
	}

	string bidAssNames = "ABCDEFGHIJKL";
	string bidSSNames = "ABCDEF";
	string monoNames = "abcdefghijkl";
	vector<string> atomsCodesLetters(numberCodes.size());
	int k = 0;
	int oldNumberCode = 0;
	int firstBid, secondBid;
	bool haveBidentate = false;
	vector<int> bidAssymetricOrdered;
	vector<int> bidAssymetricType;
	for (size_t i = 0; i < bidentateChosen.size(); i += 2)
	{
		int bidI1 = bidentateChosen[i];
		int bidI2 = bidentateChosen[i + 1];
		if (numberCodes[bidI1] != numberCodes[bidI2])
		{
			haveBidentate = true;
			if(numberCodes[bidI1] < numberCodes[bidI2])
			{
				firstBid = bidI1;
				secondBid = bidI2;
			}
			else
			{
				firstBid = bidI2;
				secondBid = bidI1;
			}
			bidAssymetricOrdered.push_back(firstBid);
			bidAssymetricOrdered.push_back(secondBid);
			if (oldNumberCode != numberCodes[firstBid])
			{
				k += 2;
				oldNumberCode = numberCodes[firstBid];
			}
			bidAssymetricType.push_back(k);
//			atomsCodesLetters[firstBid] = bidAssNames[k];
//			atomsCodesLetters[secondBid] = bidAssNames[k+1];
		}
	}
	if (haveBidentate)
	{
		oldNumberCode += 2;
		int assOld = bidAssymetricType[bidAssymetricType.size() - 1];
		int assCount = 0;
		int bidPosCount = 0;
		for (int i = (int)bidAssymetricType.size() - 1; i >= 0; i--)
		{
			if (bidAssymetricType[i] != assOld)
			{
				assOld = bidAssymetricType[i];
				assCount += 2;
			}
			int bidPos1 = bidAssymetricOrdered[bidPosCount];
			int bidPos2 = bidAssymetricOrdered[bidPosCount + 1];
			bidPosCount += 2;
			atomsCodesLetters[bidPos1] = bidAssNames[assCount];
			atomsCodesLetters[bidPos2] = bidAssNames[assCount + 1];
		}
	}

	k = 0;
	haveBidentate = false;
	vector<int> bidSymmetricType, bidSymmetricOrdered;
	for (size_t i = 0; i < bidentateChosen.size(); i += 2)
	{
		int bidI1 = bidentateChosen[i];
		int bidI2 = bidentateChosen[i + 1];
		if (numberCodes[bidI1] == numberCodes[bidI2])
		{
			haveBidentate = true;
			if (numberCodes[bidI1] != oldNumberCode)
			{
				k++;
				oldNumberCode = numberCodes[bidI1];
			}
			bidSymmetricType.push_back(k);
			bidSymmetricOrdered.push_back(bidI1);
			bidSymmetricOrdered.push_back(bidI2);
//			atomsCodesLetters[bidI1] = bidSSNames[k];
//			atomsCodesLetters[bidI2] = bidSSNames[k];
		}
	}
	if (haveBidentate)
	{
		oldNumberCode++;
		int symOld = bidSymmetricType[bidSymmetricType.size() - 1];
		int symCount = 0;
		int bidPosCount = 0;
		for (int i = (int)bidSymmetricType.size() - 1; i >= 0; i--)
		{
			if (bidSymmetricType[i] != symOld)
			{
				symOld = bidSymmetricType[i];
				symCount++;
			}
			int bidPos1 = bidSymmetricOrdered[bidPosCount];
			int bidPos2 = bidSymmetricOrdered[bidPosCount + 1];
			bidPosCount += 2;
//			atomsCodesLetters[bidPos1] = bidSSNames[symCount];
//			atomsCodesLetters[bidPos2] = bidSSNames[symCount];
		}
	}
	vector<int> monoTypes, monoPos;
	for (size_t i = 0; i < numberCodes.size(); i++)
	{
		if (numberCodes[i] >= oldNumberCode)
		{
			monoTypes.push_back(numberCodes[i] - oldNumberCode);
			monoPos.push_back(i);
//			atomsCodesLetters[i] = monoNames[numberCodes[i] - oldNumberCode];
		}
	}
	if (monoTypes.size() != 0) // ordenar esse cara junto com as posicoes
	{
		AuxMath auxMath_;
		vector<int> instructions = auxMath_.vector_ordering(monoTypes);
		auxMath_.vector_ordering_with_instructions(monoPos, instructions);
		int monoOld = monoTypes[0];
		int monoCount = 0;
		int monoPosCount = 0;
		for (int i = 0; i < (int)monoTypes.size(); i++)
		{
			if (monoTypes[i] != monoOld)
			{
				monoOld = monoTypes[i];
				monoCount++;
			}
			int monoPosI = monoPos[monoPosCount];
			monoPosCount++;
			atomsCodesLetters[monoPosI] = monoNames[monoCount];
		}
	}

	stringstream finalCode;
	for (size_t i = 0; i < atomsCodesLetters.size(); i++)
		finalCode << atomsCodesLetters[i] << " ";
	if (bidentateChosen.size() > 0)
	{
		finalCode << " / ";
		for (size_t i = 0; i < bidentateChosen.size(); i++)
			finalCode << (bidentateChosen[i] + 1) << " ";
	}
	atomTypes = numberCodes;
	for (size_t i = 0; i < atomTypes.size(); i++)
		atomTypes[i]++;
	for (size_t i = 0; i < bidentateChosen.size(); i++)
		bidentateChosen[i]++;

	return finalCode.str();
}

void ReadWriteFormats::typeLineToNumberCodes(
	string typeLines,
	std::vector<int> & atomTypes,
	std::vector<int> & bidentateChosen)
{
	string elemFirst = "ACEGIK";
	string monoStrings = "abcdefghijkl";
	string atomsCodes;
	stringstream bidLine;
	int bidNumber = 0;
	for (size_t i = 0; i < typeLines.size(); i++)
	{
		if (typeLines[i] == '/')
		{
			for (size_t j = i + 1; j < typeLines.size(); j++)
			{
				if (typeLines[j] != ' ')
				{
					bidLine << typeLines[j] << " ";
					bidNumber++;
				}
			}
			break;
		}
		if (typeLines[i] != ' ')
			atomsCodes.push_back(typeLines[i]);
	}
	for (int i = 0; i < bidNumber; i++)
	{
		int auxBid;
		bidLine >> auxBid;
		auxBid--;
		bidentateChosen.push_back(auxBid);
	}
	int k = 1;
	bool haveBidentate = false;
	vector<int> atomsCodesNumbers(atomsCodes.size());
	char oldCode = 'A';
	for (size_t i = 0; i < bidentateChosen.size(); i += 2)
	{
		int bidI1 = bidentateChosen[i];
		int bidI2 = bidentateChosen[i + 1];
		if (atomsCodes[bidI1] != atomsCodes[bidI2])
		{
			haveBidentate = true;
			size_t found = elemFirst.find(atomsCodes[bidI1]);
			int firstBid, secondBid;
			if (found != std::string::npos)
			{
				firstBid = bidI1;
				secondBid = bidI2;
			}
			else
			{
				firstBid = bidI2;
				secondBid = bidI1;
			}
			if (oldCode != atomsCodes[firstBid])
			{
				k += 2;
				oldCode = atomsCodes[firstBid];
			}
			atomsCodesNumbers[firstBid] = k;
			atomsCodesNumbers[secondBid] = k + 1;
		}
	}
	oldCode = 'A';
	if (haveBidentate)
		k += 2;
	haveBidentate = false;
	for (size_t i = 0; i < bidentateChosen.size(); i += 2)
	{
		int bidI1 = bidentateChosen[i];
		int bidI2 = bidentateChosen[i + 1];
		if (atomsCodes[bidI1] == atomsCodes[bidI2])
		{
			haveBidentate = true;
			if (!(atomsCodes[bidI1] == oldCode))
			{
				k++;
				oldCode = atomsCodes[bidI1];
			}
			atomsCodesNumbers[bidI1] = k;
			atomsCodesNumbers[bidI2] = k;
		}
	}
	if (haveBidentate)
		k++;
	for (size_t i = 0; i < atomsCodes.size(); i++)
	{
		size_t found = monoStrings.find(atomsCodes[i]);
		if (found != std::string::npos)
		{
			atomsCodesNumbers[i] = k + found;
		}
	}
	atomTypes = atomsCodesNumbers;
	for (size_t i = 0; i < bidentateChosen.size(); i++)
		bidentateChosen[i]++;
}


void ReadWriteFormats::takeRcwVgroupPointGroup(
	string line,
	int &rcw,
	string & vGroup,
	string & pGroup)
{
	stringstream convert;
	convert << line;
	convert >> rcw;
	size_t sep1Temp = line.find(";");
	size_t sep2Temp = line.find(";", sep1Temp + 1, 1);
	size_t sepGroup1 = line.find("]");
	size_t sepGroup2 = line.find("[", sepGroup1 + 1, 1);
	vGroup = line.substr(sep1Temp + 2, sep2Temp - sep1Temp - 3);
	pGroup = line.substr(sepGroup1 + 2, sepGroup2 - sepGroup1 - 3);
}

void ReadWriteFormats::takeRcwVgroupPointGroupNca(
	string line,
	int & rcw,
	int & chiral,
	int & achiral,
	string & vGroup,
	string & pGroup)
{
	stringstream convert;
	convert << line;
	convert >> rcw;
	size_t sep1Temp = line.find(";");
	size_t sep2Temp = line.find(";", sep1Temp + 1, 1);
	size_t sepGroup1 = line.find("]");
	size_t sepGroup2 = line.find("[", sepGroup1 + 1, 1);
	vGroup = line.substr(sep1Temp + 2, sep2Temp - sep1Temp - 3);
	pGroup = line.substr(sepGroup1 + 2, sepGroup2 - sepGroup1 - 3);
	if (vGroup[vGroup.size() - 2] == 'c')
	{
		chiral = 1;
		achiral = 0;
	}
	else if (vGroup[vGroup.size() - 2] == 'a')
	{
		chiral = 0;
		achiral = 1;
	}
}


void ReadWriteFormats::takeAllElementsFromCode(
	string line,
	int & rcw,
	int & chiral,
	int & achiral,
	string & vGroup,
	string & pGroup,
	string & permut)
{
	stringstream convert;
	convert << line;
	convert >> rcw;
	size_t sep1Temp = line.find(";");
	size_t sep2Temp = line.find(";", sep1Temp + 1, 1);
	size_t sepGroup1 = line.find("]");
	size_t sepGroup2 = line.find("[", sepGroup1 + 1, 1);
	size_t sepGroup3 = line.find("]", sepGroup2 + 1, 1);
	vGroup = line.substr(sep1Temp + 2, sep2Temp - sep1Temp - 3);
	pGroup = line.substr(sepGroup1 + 2, sepGroup2 - sepGroup1 - 3);
	permut = line.substr(sepGroup2 + 1, sepGroup3 - sepGroup2 - 1);


	if (vGroup[vGroup.size() - 2] == 'c')
	{
		chiral = 1;
		achiral = 0;
	}
	else if (vGroup[vGroup.size() - 2] == 'a')
	{
		chiral = 0;
		achiral = 1;
	}
}



void ReadWriteFormats::takeAllElementsFromCodeNew(
	std::string line,
	int coordination,
	int & rcw,
	string & chirality,
	std::string & vGroup,
	std::string & pGroup,
	std::string & setGroup,
	std::vector<int> & permut)
{
	stringstream convert;
	convert << line;
	convert >> rcw;
	size_t sep1Temp = line.find(";");
	size_t sep2Temp = line.find(";", sep1Temp + 1, 1);
	size_t sepGroup1 = line.find("]");
	size_t sepGroup2 = line.find("[", sepGroup1 + 1, 1);
	size_t sepGroup3 = line.find("]", sepGroup2 + 1, 1);
	vGroup = line.substr(sep1Temp + 2, sep2Temp - sep1Temp - 3);
	pGroup = line.substr(sepGroup1 + 2, sepGroup2 - sepGroup1 - 3);

	stringstream convert2;
	convert2 << pGroup;
	string aux;
	convert2 >> aux >> pGroup >> chirality >> setGroup;
	string auxPermut;

	auxPermut = line.substr(sepGroup2 + 1, sepGroup3 - sepGroup2 - 1);
	stringstream convert3;
	convert3 << auxPermut;
	permut.resize(coordination);
	for (int i = 0; i < coordination; i++)
		convert3 >> permut[i];
}


void ReadWriteFormats::symmetryGroupOrdering(
	vector<int> &uniqRcw,
	vector<string> &uniqPgroup,
	vector<int> &uniqCount)
{
	bool tradePositions = true;
	int auxRcw, auxCount;
	string auxpGroup;
	while (tradePositions)
	{
		tradePositions = false;
		for (size_t i = 0; i < uniqRcw.size() - 1; i++)
		{
			bool trade = hyerarchyOrdering(
				uniqRcw[i],
				uniqRcw[i + 1],
				uniqPgroup[i],
				uniqPgroup[i + 1]);
			if (trade)
			{
				auxRcw = uniqRcw[i];
				auxpGroup = uniqPgroup[i];
				auxCount = uniqCount[i];
				uniqRcw[i] = uniqRcw[i + 1];
				uniqPgroup[i] = uniqPgroup[i + 1];
				uniqCount[i] = uniqCount[i + 1];
				uniqRcw[i + 1] = auxRcw;
				uniqPgroup[i + 1] = auxpGroup;
				uniqCount[i + 1] = auxCount;
				tradePositions = true;
			}
		}
	}
}

bool ReadWriteFormats::hyerarchyOrdering(
	int rcw1,
	int rcw2,
	string gPoint1,
	string gPoint2)
{
	return !hyerarchyOrdering(
		gPoint1,
		gPoint2);

	string pHyerarchy = "hidvs";
	if (rcw1 > rcw2)
	{
		return false;
	}
	else if (rcw1 < rcw2)
	{
		return true;
	}
	else
	{
		int k1 = pHyerarchy.size();
		int k2 = pHyerarchy.size();
		for (size_t k = 0; k < gPoint1.size(); k++)
		{
			int pos = findCharOnString(pHyerarchy, gPoint1[k]);
			if (pos < pHyerarchy.size())
			{
				k1 = pos;
				break;
			}
		}
		for (size_t k = 0; k < gPoint2.size(); k++)
		{
			int pos = findCharOnString(pHyerarchy, gPoint2[k]);
			if (pos < pHyerarchy.size())
			{
				k2 = pos;
				break;
			}
		}
		if (k1 < k2)
			return true;
		else if (k1 > k2)
			return false;
	}

	if (rcw1 == rcw2)
	{
		if ((gPoint1[0] != gPoint2[0]))
		{
			if (gPoint1[0] == 'O')
				return false;
			else if (gPoint2[0] == 'O')
				return true;
			else if (gPoint1[0] == 'T')
				return false;
			else if (gPoint2[0] == 'T')
				return true;
			else if (gPoint1[0] == 'D')
				return false;
			else if (gPoint2[0] == 'D')
				return true;
			else if (gPoint1[0] == 'S')
				return false;
			else if (gPoint2[0] == 'S')
				return true;
		}
		else if ((gPoint1[1] != 's') &
			(gPoint1[1] != 'i') &
			(gPoint2[1] != 's') &
			(gPoint2[1] != 'i'))
		{
			stringstream convert2;
			int auxN1, auxN2;
			convert2 << gPoint1[1] << " " << gPoint2[1];
			convert2 >> auxN1 >> auxN2;
			return auxN1 < auxN2;
		}
	}

	return false;
}

string ReadWriteFormats::includeGroupPoint(
	string vCode,
	string gPoint)
{
	size_t found = vCode.find("]");
	string code1 = vCode.substr(0, found + 1);
	string code2 = vCode.substr(found + 2, vCode.size() - found - 1);
	return code1 + " " + gPoint + " " + code2;
}

string ReadWriteFormats::takeGroupPoint(string vCode)
{
	size_t sepGroup1 = vCode.find("]");
	size_t sepGroup2 = vCode.find("[", sepGroup1 + 1, 1);
	string pGroup = vCode.substr(sepGroup1 + 2, sepGroup2 - sepGroup1 - 3);
	return pGroup;
}

string ReadWriteFormats::atomTypesToAtomStrings(
	vector<int> atomTypes)
{
	stringstream convert;
	for (size_t i = 0; i < atomTypes.size(); i++)
		convert << atomLabels[atomTypes[i]] << "  ";

	return convert.str();
}

vector<int> ReadWriteFormats::atomStringToAtomTypes(
	vector<string> &atomsStrings)
{
	vector<int> atomTypes;
	for (size_t i = 0; i < atomsStrings.size(); i++)
	{
		for (size_t j = 0; j < atomLabels.size(); j++)
		{
			if (atomsStrings[i] == atomLabels[j])
				atomTypes.push_back(j);
		}
	}
	return atomTypes;
}

bool ReadWriteFormats::hyerarchyOrdering(
	string gPoint1,
	string gPoint2)
{
	if ((gPoint1[0] != gPoint2[0]))
	{
		if (gPoint1[0] == 'O')
			return false;
		else if (gPoint2[0] == 'O')
			return true;
		else if (gPoint1[0] == 'T')
			return false;
		else if (gPoint2[0] == 'T')
			return true;
		else if (gPoint1[0] == 'D')
			return false;
		else if (gPoint2[0] == 'D')
			return true;
		else if (gPoint1[0] == 'S')
			return true;
		else if (gPoint2[0] == 'S')
			return false;
	}

	// mesma letra - e C ou D - verificar rotacoes
	if ((gPoint1[1] != 's') &
		(gPoint1[1] != 'i') &
		(gPoint2[1] != 's') &
		(gPoint2[1] != 'i'))
	{
		stringstream convert2;
		int auxN1, auxN2;
		convert2 << gPoint1[1] << " " << gPoint2[1];
		convert2 >> auxN1 >> auxN2;
		if (auxN1 != auxN2)
			return auxN1 < auxN2;
	}

	string pHyerarchy = "hidvs";
	int k1 = pHyerarchy.size();
	int k2 = pHyerarchy.size();
	for (size_t k = 0; k < gPoint1.size(); k++)
	{
		int pos = findCharOnString(pHyerarchy, gPoint1[k]);
		if (pos < pHyerarchy.size())
		{
			k1 = pos;
			break;
		}
	}
	for (size_t k = 0; k < gPoint2.size(); k++)
	{
		int pos = findCharOnString(pHyerarchy, gPoint2[k]);
		if (pos < pHyerarchy.size())
		{
			k2 = pos;
			break;
		}
	}
	
	return k1 > k2;
}


bool ReadWriteFormats::hyerarchyOrdering2(
	string gPoint1,
	string gPoint2)
{
	if ((gPoint1[0] != gPoint2[0]))
	{
		if (gPoint1[0] == 'O')
			return false;
		else if (gPoint2[0] == 'O')
			return true;
		else if (gPoint1[0] == 'T')
			return false;
		else if (gPoint2[0] == 'T')
			return true;
		else if (gPoint1[0] == 'D')
			return false;
		else if (gPoint2[0] == 'D')
			return true;
		else if (gPoint1[0] == 'S')
			return false;
		else if (gPoint2[0] == 'S')
			return true;
	}



	string pHyerarchy = "hidvs";
	int k1 = pHyerarchy.size();
	int k2 = pHyerarchy.size();
	for (size_t k = 0; k < gPoint1.size(); k++)
	{
		int pos = findCharOnString(pHyerarchy, gPoint1[k]);
		if (pos < pHyerarchy.size())
		{
			k1 = pos;
			break;
		}
	}
	for (size_t k = 0; k < gPoint2.size(); k++)
	{
		int pos = findCharOnString(pHyerarchy, gPoint2[k]);
		if (pos < pHyerarchy.size())
		{
			k2 = pos;
			break;
		}
	}
	if (k1 > k2)
		return true;
	else if (k1 < k2)
		return false;
	else if ((gPoint1[1] != 's') &
		(gPoint1[1] != 'i') &
		(gPoint2[1] != 's') &
		(gPoint2[1] != 'i'))
	{
		stringstream convert2;
		int auxN1, auxN2;
		convert2 << gPoint1[1] << " " << gPoint2[1];
		convert2 >> auxN1 >> auxN2;
		return auxN1 < auxN2;
	}

	return false;
}



int ReadWriteFormats::findCharOnString(
	string word,
	char refChar)
{
	int ref = word.size();
	for (size_t i = 0; i < word.size(); i++)
	{
		if (word[i] == refChar)
		{
			ref = i;
			break;
		}
	}
	return ref;



}


void ReadWriteFormats::ReplaceAll(std::string & data, char toSearch, std::string addStr)
{
	string auxData;
        for(size_t i = 0; i < data.size(); i++)
        {
                if(data[i] == toSearch)
                {
                        auxData += addStr;
                }
                auxData += data[i];
        }
	data = auxData;

}
 
bool ReadWriteFormats::isachiral(string pGroup)
{
	return ((findCharOnString(pGroup, 's')
		!= (int)pGroup.size()) ||
		(findCharOnString(pGroup, 'h')
			!= (int)pGroup.size()) ||
			(findCharOnString(pGroup, 'v')
				!= (int)pGroup.size()) ||
				(findCharOnString(pGroup, 'd')
					!= (int)pGroup.size()) ||
					(findCharOnString(pGroup, 'i')
						!= (int)pGroup.size()));
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

	// same as isomers to mol
	atomLabels.resize(12);
	atomLabels[0] = "O";// red
	atomLabels[1] = "P";// green
	atomLabels[2] = "N";// dark blue
	atomLabels[3] = "I";// purple
	atomLabels[4] = "Na";// brown
	atomLabels[5] = "B";// light blue
	atomLabels[6] = "Se";// dark green
	atomLabels[7] = "C";// grey
	atomLabels[8] = "He";// white
	atomLabels[9] = "Ca";// light green
	atomLabels[10] = "Ti";// gold
	atomLabels[11] = "Sc";// gold
	
}











