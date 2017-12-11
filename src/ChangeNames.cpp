#include "ChangeNames.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "AllMolecularFormulas.h"
#include "AuxMath.h"
#include "CauchyIndex.h"

using namespace std;

ChangeNames::ChangeNames(){}

ChangeNames::~ChangeNames(){}

void ChangeNames::changeNameOfFiles(string responseName)
{
	ifstream response_(responseName.c_str());
	string line;
	ofstream counting_("counting.csv");
	while (!response_.eof())
	{
		getline(response_, line);
		if (line == "")
			break;
		string combination;
		stringstream convert;
		convert << line;
		convert >> combination;
		AllMolecularFormulas allMol_;
		vector< vector<int> > combinationCode = allMol_.stringToNumber(combination);
		string newCombinationName = allMol_.newCodeToString(combinationCode);
		int systemSize = calculateSystemSize(combinationCode);
		string geomName = sizeToGeometryCode(systemSize);
		ofstream newFile_((geomName + "-" + newCombinationName + ".csv").c_str());
		string typeLine = generateNewTypeLine(combination, systemSize);
		newFile_ << typeLine << endl;

		ifstream isomerFile_((("final-") + combination).c_str());
		vector<vultorGroup> vGroup;
		calculateVultorGroup(isomerFile_, vGroup);

		isomerFile_.open((("final-") + combination).c_str());
		changeFormat(
			isomerFile_,
			counting_,
			newFile_,
			vGroup,
			geomName,
			newCombinationName,
			systemSize);
	}
}


vector<vultorGroup> ChangeNames::setVultorGroup(
	vector<int> &probs,
	vector<int> &number,
	vector<bool> &chiral)
{
	vector<vultorGroup> wholeGroup;
	AuxMath auxMath_;
	vector<int> instruct = auxMath_.vector_ordering(probs);
	auxMath_.vector_ordering_with_instructions(number, instruct);
	auxMath_.vector_ordering_with_instructions(chiral, instruct);
	int k = 0;
	for (size_t i = 0; i < probs.size(); i++)
	{
		vultorGroup newLetter;
		newLetter.blockName = takeLetter(k);
		k++;
		if (i + 1 == probs.size())
		{
			if (chiral[i])
			{
				newLetter.achiralN = 0;
				newLetter.achiralProb = 0;
				newLetter.chiralN = number[i];
				newLetter.chiralProb = probs[i];
			}
			else
			{
				newLetter.achiralN = number[i];
				newLetter.achiralProb = probs[i];
				newLetter.chiralN = 0;
				newLetter.chiralProb = 0;
			}
		}
		else
		{
			//check next to see if probability is equal.
			if (probs[i] == probs[i + 1])
			{
				if (chiral[i])
				{
					newLetter.achiralN = number[i + 1];
					newLetter.achiralProb = probs[i + 1];
					newLetter.chiralN = number[i];
					newLetter.chiralProb = probs[i];
					i++;
				}
				else
				{
					newLetter.achiralN = number[i];
					newLetter.achiralProb = probs[i];
					newLetter.chiralN = number[i + 1];
					newLetter.chiralProb = probs[i + 1];
					i++;
				}
			}
			else
			{

				if (chiral[i])
				{
					newLetter.achiralN = 0;
					newLetter.achiralProb = 0;
					newLetter.chiralN = number[i];
					newLetter.chiralProb = probs[i];
				}
				else
				{
					newLetter.achiralN = number[i];
					newLetter.achiralProb = probs[i];
					newLetter.chiralN = 0;
					newLetter.chiralProb = 0;
				}
			}
		}
		wholeGroup.push_back(newLetter);
	}

	return wholeGroup;
}

string ChangeNames::takeLetter(int nGroup)
{
	switch (nGroup)
	{
	case 0:
		return "A";
	case 1:
		return "B";
	case 2:
		return "C";
	case 3:
		return "D";
	case 4:
		return "E";
	case 5:
		return "F";
	case 6:
		return "G";
	case 7:
		return "H";
	case 8:
		return "I";
	case 9:
		return "J";
	case 10:
		return "K";
	case 11:
		return "L";
	case 12:
		return "M";
	case 13:
		return "N";
	case 14:
		return "O";
	case 15:
		return "P";
	case 16:
		return "Q";
	case 17:
		return "R";
	case 18:
		return "S";
	case 19:
		return "T";
	case 20:
		return "U";
	case 21:
		return "V";
	case 22:
		return "W";
	case 23:
		return "X";
	case 24:
		return "Y";
	case 25:
		return "Z";
	default:
		cout << "ERROR ON takeLetter" << endl;
		exit(1);
		break;
	}
	return "";
}

vultorGroup ChangeNames::findVultorGroup(int prob, vector<vultorGroup> &group)
{
	for (size_t i = 0; i < group.size(); i++)
	{
		if (
			(group[i].achiralProb == prob) ||
			(group[i].chiralProb == prob))
		{
			return group[i];
		}
	}
	cout << "findVultorGroup not found" << endl;
	exit(1);
	return group[0];
}

string ChangeNames::sizeToGeometryCode(int size)
{
	switch (size)
	{
	case 4:
		return "T-4";
		break;
	case 5:
		return "TBPY-5";
		break;
	case 6:
		return "OC-6";
		break;
	case 7:
		return "COC-7";
		break;
	case 8:
		return "SAPR-8";
		break;
	case 9:
		return "TCTPR-9";
		break;
	case 10:
		return "JMBIC-10";
		break;
	default:
		cout << "size not found" << endl;
		exit(1);
	}
}

void ChangeNames::calculateVultorGroup(
	std::ifstream & isomerFile_,
	vector<vultorGroup> &vGroup)
{
	string isomerLine;
	vector<int> vultorGroupProbs;
	vector<int> vultorGroupCounting;
	vector<bool> vultorGroupChirality; // true is chiral
	while (!isomerFile_.eof())
	{
		getline(isomerFile_, isomerLine);
		if (isomerLine == "")
			break;

		stringstream convertWeights;
		int auxWeight1, auxWeight2;
		convertWeights << isomerLine;
		convertWeights >> auxWeight1;
		auxWeight1++;
		getline(isomerFile_, isomerLine);
		if (isomerLine == "")
		{
			// not chiral
			// find value
			bool found = false;
			for (size_t i = 0; i < vultorGroupProbs.size(); i++)
			{
				if ((auxWeight1 == vultorGroupProbs[i]) && (!vultorGroupChirality[i]))
				{
					vultorGroupCounting[i]++;
					found = true;
					break;
				}
			}
			if (!found)
			{
				vultorGroupProbs.push_back(auxWeight1);
				vultorGroupCounting.push_back(1);
				vultorGroupChirality.push_back(false);
			}
		}
		else
		{
			//chiral
			stringstream convertW2;
			convertW2 << isomerLine;
			convertW2 >> auxWeight2;
			auxWeight2++;
			//chiral
			getline(isomerFile_, isomerLine);
			if (isomerLine != "")
			{
				cout << "ERROR ON CHANGE NAMES - EMPTY LINE EXPECTED" << endl;
				exit(1);
			}

			bool found = false;
			int totalWeight = auxWeight1 + auxWeight2;
			for (size_t i = 0; i < vultorGroupProbs.size(); i++)
			{
				if ((totalWeight == vultorGroupProbs[i]) && (vultorGroupChirality[i]))
				{
					vultorGroupCounting[i]++;
					found = true;
					break;
				}
			}
			if (!found)
			{
				vultorGroupProbs.push_back(totalWeight);
				vultorGroupCounting.push_back(1);
				vultorGroupChirality.push_back(true);
			}
		}
	}
	isomerFile_.close();

	vector<vultorGroup> auxVGroup = setVultorGroup(
		vultorGroupProbs,
		vultorGroupCounting,
		vultorGroupChirality);
	vGroup = auxVGroup;

	vector<double> entropyOrdering(auxVGroup.size());
	for (size_t i = 0; i < auxVGroup.size(); i++)
	{
		entropyOrdering[i] = auxVGroup[i].achiralN * auxVGroup[i].achiralProb
			+ 2 * auxVGroup[i].chiralN * auxVGroup[i].chiralProb;
	}
	AuxMath auxMath_;
	vector<int> instructions = auxMath_.vector_ordering(entropyOrdering);
	// ORDER VECTOR WITH INSTRUCTIONS
	for (int i = 0; i < instructions.size(); i += 2)
	{
		cout << "inversion" << endl;
		vGroup[instructions[i]] = auxVGroup[instructions[i + 1]];
		vGroup[instructions[i + 1]] = auxVGroup[instructions[i]];
	}
}


void ChangeNames::changeFormat(
	ifstream &isomerFile_,
	ofstream & counting_,
	ofstream & newFile_,
	vector<vultorGroup> &vGroup,
	string & geomName,
	string & newCombinationName,
	int systemSize
)
{
	string isomerLine;
	int totalChiral = 0;
	int totalAchiral = 0;
	bool chiral;
	CauchyIndex cauchy_(systemSize);
	while (!isomerFile_.eof())
	{
		getline(isomerFile_, isomerLine);
		if (isomerLine == "")
		{
			newFile_ << endl;
			continue;
		}

		//check next to see if chiral
		string isomerLine2;
		if (isomerFile_.eof())
			chiral = false;
		else
		{
			getline(isomerFile_, isomerLine2);
			if (isomerLine2 == "")
				chiral = false;
			else
				chiral = true;
		}
		stringstream convertLine1;
		int auxWeight;
		string dummy1, dummy2, dummy3;
		convertLine1 << isomerLine;
		convertLine1 >> auxWeight;
		convertLine1 >> dummy1 >> dummy2 >> dummy3;
		vector<int> notation1(systemSize);
		for (int i = 0; i < systemSize; i++)
			convertLine1 >> notation1[i];

		newFile_ << (auxWeight + 1) << " ; "
			<< "{" << geomName << " "
			<< "[" << (notation1[0] + 1);
		for (size_t i = 1; i < notation1.size(); i++)
			newFile_ << " " << (notation1[i] + 1);
		newFile_ << "] ";


		if (!chiral)
		{
			vultorGroup thisGroup = findVultorGroup(auxWeight + 1, vGroup);
			newFile_ << thisGroup.blockName << "a}" << endl << endl;
			totalAchiral++;
		}
		else
		{
			vultorGroup thisGroup = findVultorGroup(2 * (auxWeight + 1), vGroup);
			totalChiral++;
			newFile_ << thisGroup.blockName << "c}" << endl;
			stringstream convertLine2;
			int auxWeight2;
			convertLine2 << isomerLine2;
			convertLine2 >> auxWeight2;
			convertLine2 >> dummy1 >> dummy2 >> dummy3;
			vector<int> notation2(systemSize);
			for (int i = 0; i < systemSize; i++)
				convertLine2 >> notation2[i];

			newFile_ << (auxWeight2 + 1) << " ; "
				<< "{" << geomName << " "
				<< "[" << (notation2[0] + 1);
			for (size_t i = 1; i < notation2.size(); i++)
				newFile_ << " " << (notation2[i] + 1);
			newFile_ << "] ";
			newFile_ << thisGroup.blockName << "c}" << endl;
		}
	}
	counting_ << newCombinationName << " ; "
		<< totalChiral << " ; "
		<< totalAchiral << " ; ";
	for (size_t i = 0; i < vGroup.size(); i++)
	{
		if (vGroup[i].chiralN != 0)
		{
			counting_ << vGroup[i].blockName << "c ; "
				<< vGroup[i].chiralN << " ; "
				<< vGroup[i].chiralProb << " ; ";
		}
		if (vGroup[i].achiralN != 0)
		{
			counting_ << vGroup[i].blockName << "a ; "
				<< vGroup[i].achiralN << " ; "
				<< vGroup[i].achiralProb << " ; ";
		}
	}
	counting_ << endl;
}


int ChangeNames::calculateSystemSize(std::vector< std::vector<int> > & combinationCode)
{
	int systemSize = 0;
	for (size_t i = 0; i < combinationCode.size(); i++)
	{
		for (size_t j = 0; j < combinationCode[i].size(); j++)
		{
			if (i > 0)
				systemSize += 2 * combinationCode[i][j];
			else
				systemSize += combinationCode[i][j];
		}
	}
	return systemSize;
}


string ChangeNames::generateNewTypeLine(string &combination, int systemSize)
{
	ifstream typesFile_((combination + "---atomTypes.txt").c_str());
	string typeLine;
	getline(typesFile_, typeLine);
	stringstream convert;
	convert << typeLine;
	vector<int> composition(systemSize);
	for (int i = 0; i < systemSize; i++)
	{
		convert >> composition[i];
	}
	string dummy;
	convert >> dummy;
	convert >> dummy;
	vector<int> bidentates;
	int readingBid;
	while (true)
	{
		convert >> readingBid;
		if (readingBid < 0)
			break;
		else
			bidentates.push_back(readingBid);
	}
	stringstream convert2;
	for (int i = 0; i < systemSize; i++)
	{
		convert2 << composition[i] + 1  << "  ";
	}
	if (bidentates.size() != 0)
	{
		convert2 << "/  ";
		for (int i = 0; i < bidentates.size(); i++)
		{
			convert2 << bidentates[i] + 1 << "  ";
		}
	}
	return convert2.str();
}
