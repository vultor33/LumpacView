#include "ChangeNames.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <algorithm>
#include <iterator>

#include "AllMolecularFormulas.h"
#include "AuxMath.h"
#include "Geometries.h"
#include "ReadWriteFormats.h"

using namespace std;

ChangeNames::ChangeNames(){}

ChangeNames::~ChangeNames(){}

void ChangeNames::changeNameOfFiles(
	string responseName,
	int geoCode,
	string pathRead,
	string pathWrite)
{
	Geometries geo_;
	string geomName = geo_.sizeToGeometryCode(geoCode);
	ifstream response_((pathRead + responseName).c_str());
	string line;
	string countingName = pathWrite + geomName + "-counting.csv";
	ofstream counting_(countingName.c_str());

	int systemSize = 0;
	vector<string> allTypeLinesTemp;
	vector<string> oldCombinationNamesTemp;
	vector<string> allNewCombinationNamesTemp;
	vector < vector< vector<int> > > allCombCodes;
	while (!response_.eof())
	{
		getline(response_, line);
		if (line == "")
			break;
		string combination;
		stringstream convert;
		convert << line;
		convert >> combination;
		ReadWriteFormats rwf_;
		vector< vector<int> > combinationCode = rwf_.compositionToNumberOld(combination);
		string newCombinationName = rwf_.newCodeToString(combinationCode);
		newCombinationName = "M" + newCombinationName;
		if (systemSize == 0)
			systemSize = calculateSystemSize(combinationCode);
		string typeLine = generateNewTypeLine(pathRead, combination, systemSize);
		allTypeLinesTemp.push_back(typeLine);
		oldCombinationNamesTemp.push_back(combination);
		allNewCombinationNamesTemp.push_back(newCombinationName);
		allCombCodes.push_back(combinationCode);
	}

	//order combCodes
	size_t totalAcceptSize = 0;
	vector<string> allTypeLines;
	vector<string> oldCombinationNames;
	vector<string> allNewCombinationNames;
	for (size_t k = 0; k < 6; k++)
	{
		for (int i = 0; i < (int)allCombCodes.size(); i++)
		{
			if ((allCombCodes[i][1].size() + allCombCodes[i][2].size()) == k)
			{
				allTypeLines.push_back(allTypeLinesTemp[i]);
				oldCombinationNames.push_back(oldCombinationNamesTemp[i]);
				allNewCombinationNames.push_back(allNewCombinationNamesTemp[i]);
			}
		}
	}

	// allTypeLines
	// tem bidentados?
	// os bidentados sao iguais ou diferentes?
	// se forem diferentes,
	// pegue o menor numero, coloque A / B / C / D e etc.
	// ande em tudo, o mesmo numero e o mesmo código.
	// ande no codigo dos assimétricos.
	// se forem iguais acabou o assimétricos. comece o codigo dos simetricos.
	// e assim por diante.

	// allTypeLines reverse.
	// 
	/* E POSSIVEL FAZER A REVERSAO MAS ELE MUDA A COMPOSICAO DA PARADA
	ReadWriteFormats rwf_;
	for (int i = 0; i < allTypeLines.size(); i++)
	{
		vector<int> types1, types2, bid1, bid2;
		string newName = rwf_.typeLineToLetters(
			allTypeLines[i],
			types1,
			bid1);
		rwf_.typeLineToNumberCodes(
			newName,
			types2,
			bid2);
		if ((types1 != types2) || (bid1 != bid2))
		{
			cout << "ChangeNames::changeNameOfFiles crytographic error" << endl;
			exit(1);
		}
		allTypeLines[i] = newName;
	}
	*/


	for(size_t i = 0; i < oldCombinationNames.size(); i++)
	{
		ifstream isomerFile_((pathRead + "final-" + oldCombinationNames[i]).c_str());
		vector<vultorGroup> vGroup;
		calculateVultorGroup(isomerFile_, vGroup);
		isomerFile_.close();
		isomerFile_.open((pathRead + "final-" + oldCombinationNames[i]).c_str());
		string geomCombName = geomName + " [" + allNewCombinationNames[i] + "]";
		ofstream newFile_((pathWrite + geomName + "-" + allNewCombinationNames[i] + ".csv").c_str());
		newFile_ << allTypeLines[i] << endl;
		cout << "i:  " << i << "  " << geomCombName << endl;
		changeFormat(
			isomerFile_,
			counting_,
			newFile_,
			vGroup,
			geomCombName,
			allNewCombinationNames[i],
			systemSize);
	}
}


void ChangeNames::createNewCounting(
	int geoCode,
	string pathRead,
	string responseName)
{
	Geometries geo_;
	string geomName = geo_.sizeToGeometryCode(geoCode);

	ReadWriteFormats rwf_;
	ifstream response_((pathRead + responseName).c_str());
	string line;

	vector<string> allFormulas;
	vector< vector<int> > allRcw;
	vector< vector<int> > allCount;
	vector< vector<string> > allPgroup;
	vector< vector<int> > allRce;
	vector< vector<int> > allNumbersCA;
	vector< vector<string> > allLettersSets;

	while (!response_.eof())
	{
		getline(response_, line);
		if (line == "")
			break;
		string combination;
		stringstream convert;
		convert << line;
		convert >> combination;
		vector< vector<int> > combinationCode = rwf_.compositionToNumberOld(combination);
		string newCombinationName = rwf_.newCodeToString(combinationCode);
		newCombinationName = "M" + newCombinationName;
		string allIsomersCombinationFile = geomName + "-" + newCombinationName + ".csv";

		vector<int> rcw;
		vector<int> count;
		vector<string> pGroup;
		vector<int> numbersCA;
		string firstLine;
		vector<string> listOfPermutations;
		vector<string> listOfChiralities;
		std::vector<int> listOfRcw;
		std::vector<std::string> listOfPGroup;

		generateOrderingGroupPoint(
			pathRead + allIsomersCombinationFile,
			rcw,
			count,
			numbersCA,
			pGroup,
			firstLine,
			listOfPermutations,
			listOfChiralities,
			listOfRcw,
			listOfPGroup);
			

	/*	
		generateOrderingGroupPoint(
			pathRead + allIsomersCombinationFile,
			rcw,
			count,
			numbersCA,
			pGroup);
	*/

		vector<int> rce(rcw.size());
		for (size_t i = 0; i < rcw.size(); i++)
			rce[i] = rcw[i] * count[i];

		vector<string> lettersSetsGroups;
		rceOrderingCriteria(rcw, rce, count, pGroup, lettersSetsGroups);

		allRcw.push_back(rcw);
		allCount.push_back(count);
		allNumbersCA.push_back(numbersCA);
		allPgroup.push_back(pGroup);
		allFormulas.push_back(newCombinationName);
		allRce.push_back(rce);
		allLettersSets.push_back(lettersSetsGroups);

		ofstream printIso_(("resultados/" + allIsomersCombinationFile).c_str());
		printIso_ << firstLine << endl;
		for (size_t i = 0; i < listOfPermutations.size(); i++)
		{
			string letterEachIsomer = "";
			for (size_t j = 0; j < pGroup.size(); j++)
			{
				if (pGroup[j] == listOfPGroup[i])
				{
					letterEachIsomer = lettersSetsGroups[j];
					break;
				}
			}
			if (letterEachIsomer == "")
			{
				cout << "ERROR ON: ChangeNames::createNewCounting" << endl;
				exit(1);
			}

			//puxar informacao pela simetria
			/*
			listOfPermutations,
				listOfChiralities,
				listOfRcw,
				listOfPGroup);
				geomName
				newCombinationName
				lettersSetsGroups
				*/

			printIso_ << listOfRcw[i]
				<< " ; {["
				<< newCombinationName
				<< "] "
				<< geomName
				<< " "
				<< listOfPGroup[i]
				<< " "
				<< listOfChiralities[i]
				<< " "
				<< letterEachIsomer
				<< " ["
				<< listOfPermutations[i]
				<< "]}"
				<< endl;
			if (listOfChiralities[i] == "c")
			{
				i++;
				printIso_ << listOfRcw[i]
					<< " ; {["
					<< newCombinationName
					<< "] "
					<< geomName
					<< " "
					<< listOfPGroup[i]
					<< " "
					<< listOfChiralities[i]
					<< " "
					<< letterEachIsomer
					<< " ["
					<< listOfPermutations[i]
					<< "]}"
					<< endl;
			}
			printIso_ << endl;
		}
		printIso_.close();

	}



	/* reordenar com hierarquia
	bool tradePositions = true;
	vector<int> auxRcw;
	vector<int> auxCount;
	vector<string> auxPgroup;
	string auxFormula;
	while (tradePositions)
	{
		tradePositions = false;
		for (size_t i = 0; i < allRcw.size() - 1; i++)
		{
			bool trade = false;
			if (allRcw[i].size() > allRcw[i + 1].size())
				continue;
			else if ((allRcw[i].size() < allRcw[i + 1].size()))
				trade = true;
			else
			{
				for (size_t k = 0; k < allRcw[i].size(); k++)
				{
					if (allPgroup[i][k] != allPgroup[i + 1][k])
					{
						bool hyearc = rwf_.hyerarchyOrdering(
							allPgroup[i][k],
							allPgroup[i + 1][k]);
						if (!hyearc)
						{
							trade = true;
							break;
						}
						else
							break;
					}
				}
			}
			if (trade)
			{
				auxRcw = allRcw[i];
				auxPgroup = allPgroup[i];
				auxCount = allCount[i];
				auxFormula = allFormulas[i];
				allRcw[i] = allRcw[i + 1];
				allPgroup[i] = allPgroup[i + 1];
				allCount[i] = allCount[i + 1];
				allFormulas[i] = allFormulas[i + 1];
				allRcw[i + 1] = auxRcw;
				allPgroup[i + 1] = auxPgroup;
				allCount[i + 1] = auxCount;
				allFormulas[i + 1] = auxFormula;
				tradePositions = true;
			}
		}
	}
	*/

	ofstream counting_((geomName + "-counting.csv").c_str());

	// Complete header
	vector<string> header = allPgroup[0];
	for (size_t i = 1; i < allPgroup.size(); i++)
	{
		for (size_t j = 0; j < allPgroup[i].size(); j++)
		{
			vector<string>::iterator it = find(header.begin(), header.end(), allPgroup[i][j]);
			if (it == header.end())
			{
				header.push_back(allPgroup[i][j]);
			}
		}
	}

	/*
	tradePositions = true;
	while (tradePositions)
	{
		string auxHeaderTrade;
		tradePositions = false;
		for (size_t i = 0; i < header.size() - 1; i++)
		{
			bool hyearc = rwf_.hyerarchyOrdering(
				header[i],
				header[i+1]);
			if (!hyearc)
			{
				auxHeaderTrade = header[i];
				header[i] = header[i + 1];
				header[i + 1] = auxHeaderTrade;
				tradePositions = true;
			}
		}
	}
	*/

	/* HEADER 
	vector<string> auxPgroup = vector<string>();
	//reverse(header.begin(), header.end());
	counting_ << "Formula;RCWP;";
	for (size_t i = 0; i < header.size(); i++)
		counting_ << header[i] << ";";
	counting_ << endl;
	for (size_t i = 0; i < allPgroup.size(); i++)
	{
		vector<int> rce(allRcw[i].size());
		for (size_t l = 0; l < allRcw[i].size(); l++)
			rce[l] = allCount[i][l] * allRcw[i][l];
		int rceMin = *min_element(rce.begin(), rce.end());
		int rcwMin = *min_element(allRcw[i].begin(), allRcw[i].end());

		counting_ << allFormulas[i] << "; ";

		// RCP COMENTAR
		for (size_t k = 0; k < allRcw[i].size() - 1; k++)
		{
			if (rce[k] % rceMin != 0)
				counting_ << fixed << setprecision(1) << (double)rce[k] / (double)rceMin << ":";
			else
				counting_ << rce[k] / rceMin << ":";
		}
		if (rce[rce.size() - 1] % rceMin != 0)
			counting_ << fixed << setprecision(1) << (double)rce[rce.size() - 1] / (double)rceMin;
		else
			counting_ << rce[rce.size() - 1] / rceMin;
		counting_ << "; ";
		//////////////////////////////////////////////

		// RCW
		for (size_t k = 0; k < allRcw[i].size() - 1; k++)
		{
			if (allRcw[i][k] % rcwMin != 0)
				counting_ << fixed << setprecision(1) << (double)allRcw[i][k] / (double)rcwMin << ":";
			else
				counting_ << allRcw[i][k] / rcwMin << ":";
		}
		if (allRcw[i][allRcw[i].size() - 1] % rcwMin != 0)
			counting_ << fixed << setprecision(1) << (double)allRcw[i][allRcw[i].size() - 1] / (double)rcwMin;
		else
			counting_ << allRcw[i][allRcw[i].size() - 1] / rcwMin;
		counting_ << ";";


		for (size_t j = 0; j < header.size(); j++)
		{
			bool marked = false;
			for (size_t k = 0; k < allPgroup[i].size(); k++)
			{
				if (allPgroup[i][k] == header[j])
				{
					counting_ << allCount[i][k] << ";";
					marked = true;
					break;
				}
			}
			if (!marked)
				counting_ << "-;";
		}
		counting_ << endl;
	}
	*/


	/* Own Symmetry groups 	
	//ofstream cheackOrder_("checkOrder.csv");
	for (size_t i = 0; i < allPgroup.size(); i++)
	{
		//for (size_t k = 0; k < allPgroup[i].size(); k++)
		//{
		//	cheackOrder_ << allPgroup[i][k] << " ; ";
		//}
		//cheackOrder_ << endl;

		if (allPgroup[i] != auxPgroup)
		{
			counting_ << endl;
			auxPgroup = allPgroup[i];
			counting_ << "Formula;RCP;RCW;";
			for (size_t k = 0; k < allPgroup[i].size(); k++)
				counting_ << allPgroup[i][k] << ";";
			counting_ << endl;
		}

		counting_ << allFormulas[i] << "; ";
		vector<int> rce(allRcw[i].size());
		for (size_t l = 0; l < allRcw[i].size(); l++)
			rce[l] = allCount[i][l] * allRcw[i][l];
		int rceMin = *min_element(rce.begin(), rce.end());
		int rcwMin = *min_element(allRcw[i].begin(), allRcw[i].end());
		// RCP
		for (size_t k = 0; k < allRcw[i].size() - 1; k++)
		{
			if (rce[k] % rceMin != 0)
				counting_ << fixed << setprecision(1) << (double)rce[k] / (double)rceMin << ":";
			else
				counting_ << rce[k] / rceMin << ":";
		}
		if (rce[rce.size() - 1] % rceMin != 0)
			counting_ << fixed << setprecision(1) << (double)rce[rce.size() - 1] / (double)rceMin;
		else
			counting_ << rce[rce.size() - 1] / rceMin;
		counting_ << "; ";


		// RCW
		for (size_t k = 0; k < allRcw[i].size() - 1; k++)
		{
			if (allRcw[i][k] % rcwMin != 0)
				counting_ << fixed << setprecision(1) << (double)allRcw[i][k] / (double)rcwMin << ":";
			else
				counting_ << allRcw[i][k] / rcwMin << ":";
		}
		if (allRcw[i][allRcw[i].size() - 1] % rcwMin != 0)
			counting_ << fixed << setprecision(1) << (double)allRcw[i][allRcw[i].size() - 1] / (double)rcwMin;
		else
			counting_ << allRcw[i][allRcw[i].size() - 1] / rcwMin;
		counting_ << ";";

		for (size_t k = 0; k < allCount[i].size(); k++)
			counting_ << allCount[i][k] << ";";

		counting_ << endl;
	}
	*/




/*
	for (size_t i = 0; i < allPgroup.size(); i++)
	{
		vector<int> rce(allRcw[i].size());
		for (size_t l = 0; l < allRcw[i].size(); l++)
			rce[l] = allCount[i][l] * allRcw[i][l];
		int rceMin = *min_element(rce.begin(), rce.end());
		int rcwMin = *min_element(allRcw[i].begin(), allRcw[i].end());
		counting_ << allFormulas[i] << "; ";

		// RCP
		for (size_t k = 0; k < allRcw[i].size() - 1; k++)
		{
			if (rce[k] % rceMin != 0)
				counting_ << fixed << setprecision(1) << (double)rce[k] / (double)rceMin << ":";
			else
				counting_ << rce[k] / rceMin << ":";
		}
		if (rce[rce.size() - 1] % rceMin != 0)
			counting_ << fixed << setprecision(1) << (double)rce[rce.size() - 1] / (double)rceMin;
		else
			counting_ << rce[rce.size() - 1] / rceMin;
		counting_ << "; ";

		// RCWP
		for (size_t k = 0; k < allRcw[i].size() - 1; k++)
		{
			if (allRcw[i][k] % rcwMin != 0)
				counting_ << fixed << setprecision(1) << (double)allRcw[i][k] / (double)rcwMin << ":";
			else
				counting_ << allRcw[i][k] / rcwMin << ":";
		}
		if (allRcw[i][allRcw[i].size() - 1] % rcwMin != 0)
			counting_ << fixed << setprecision(1) << (double)allRcw[i][allRcw[i].size() - 1] / (double)rcwMin;
		else
			counting_ << allRcw[i][allRcw[i].size() - 1] / rcwMin;
		counting_ << ";";

		for (size_t k = 0; k < allCount[i].size(); k++)
		{
			counting_ << allPgroup[i][k] << ";"
				<< allCount[i][k] << ";";
		}

		counting_ << endl;
	}
*/


/* RCE BASED ORDERING
for (size_t i = 0; i < allPgroup.size(); i++)
{
	int rceMin = *min_element(allRce[i].begin(), allRce[i].end());
	counting_ << allFormulas[i] << "; ";

	counting_ << allNumbersCA[i][0] + allNumbersCA[i][1]
		<< ";" << allNumbersCA[i][0]
		<< ";" << allNumbersCA[i][1] << ";";

	// RCP
	for (size_t k = 0; k < allRce[i].size() - 1; k++)
	{
		if (allRce[i][k] % rceMin != 0)
		{
			counting_ << fixed 
				<< setprecision(1) 
				<< (double)allRce[i][k] / (double)rceMin
				<< "(" << allPgroup[i][k] << "):";
		}
		else
		{
			counting_ << allRce[i][k] / rceMin 
				<< "(" << allPgroup[i][k] << "):";
		}
	}
	if (allRce[i][allRce[i].size() - 1] % rceMin != 0)
	{
		counting_ << fixed 
			<< setprecision(1) 
			<< (double)allRce[i][allRce[i].size() - 1] / (double)rceMin
			<< "(" << allPgroup[i][allRce[i].size() - 1] << ")";
	}
	else
	{
		counting_ << allRce[i][allRce[i].size() - 1] / rceMin
			<< "(" << allPgroup[i][allRce[i].size() - 1] << ")";
	}
	counting_ << endl;
}
*/

for (size_t i = 0; i < allPgroup.size(); i++)
{
	int rceMin = *min_element(allRce[i].begin(), allRce[i].end());
	counting_ << allFormulas[i] << ";";

	counting_ << allNumbersCA[i][0] + allNumbersCA[i][1]
		<< ";" << allNumbersCA[i][0]
		<< ";" << allNumbersCA[i][1] << ";";
	
	counting_ << " ";

	// RCP
	for (size_t k = 0; k < allRce[i].size() - 1; k++)
	{
		if (allRce[i][k] % rceMin != 0)
		{
			counting_ << fixed
				<< setprecision(1)
				<< (double)allRce[i][k] / (double)rceMin
				<< ":";
		}
		else
		{
			counting_ << allRce[i][k] / rceMin
				<< ":";
		}
	}
	if (allRce[i][allRce[i].size() - 1] % rceMin != 0)
	{
		counting_ << fixed
			<< setprecision(1)
			<< (double)allRce[i][allRce[i].size() - 1] / (double)rceMin;
	}
	else
	{
		counting_ << allRce[i][allRce[i].size() - 1] / rceMin;
	}
	counting_ << ";";

	int k = 0;
	for (size_t j = 0; j < allRce[i].size(); j++)
	{
		k++;
		if (k == 3)
		{
			counting_ << endl;
			counting_ << ";;;;;";
			k = 0;
		}
		
		string chiralLetter;
		if (rwf_.isachiral(allPgroup[i][j]))
			chiralLetter = "a";
		else
			chiralLetter = "c";

		counting_ << allLettersSets[i][j] << ";"
			<< allPgroup[i][j] << ";"
			<< chiralLetter << ";"
			<< allCount[i][j] << ";"
			<< allRcw[i][j] << ";";
	}
	counting_ << endl;
}



}



void ChangeNames::rceOrderingCriteria(
	std::vector<int> &rcw,
	std::vector<int> &rce,
	std::vector<int> &count,
	std::vector<std::string> &pGroup,
	std::vector<std::string> &setsLetter)
{
	// Transformar em uma function
	// REORDENAR AQUI COM RCP E CRITERIOS.
	bool tradePositions;
	ReadWriteFormats rwf_;
	tradePositions = true;
	while (tradePositions)
	{
		tradePositions = false;
		for (size_t j = 0; j < (rce.size() - 1); j++)
		{
			bool trade = false;
			if (rce[j] > rce[j + 1])
				continue;
			else if (rce[j] < rce[j + 1])
				trade = true;
			else // iguais
			{
				bool hyearc = rwf_.hyerarchyOrdering(
					pGroup[j],
					pGroup[j + 1]);
				if (!hyearc)
				{
					trade = true;
				}
				else
					continue;
			}
			if (trade)
			{
				int rcwj = rcw[j];
				int rcej = rce[j];
				int contj = count[j];
				string pGroupj = pGroup[j];

				rcw[j] = rcw[j + 1];
				rce[j] = rce[j + 1];
				count[j] = count[j + 1];
				pGroup[j] = pGroup[j + 1];

				rcw[j + 1] = rcwj;
				rce[j + 1] = rcej;
				count[j + 1] = contj;
				pGroup[j + 1] = pGroupj;
				tradePositions = true;
			}
		}
	}


	vector<int> nlinhas(rce.size());
	for (size_t i = 0; i < rce.size(); i++)
		nlinhas[i] = 0;
	setsLetter.resize(rce.size());

	//linha igual repeticao da letra
	for (size_t i = 0; i < rce.size(); i++)
	{
		if (rce.size() == 1)
		{
			nlinhas[0] = 0;
			break;
		}
		// checa se o proximo e igual o anterior.
		if (i != 0)
		{
			if (rce[i - 1] == rce[i])
			{
				nlinhas[i] = nlinhas[i - 1] + 1;
			}
			else if (i != rce.size() - 1)
			{
				if (rce[i + 1] == rce[i])
				{
					nlinhas[i] = 1;
				}
			}
		}
		else
		{
			if (rce[i + 1] == rce[i])
			{
				nlinhas[i] = 1;
			}
		}
	}

	int kLetter = -1;
	setsLetter.resize(rce.size());
	for (size_t i = 0; i < rce.size(); i++)
	{
		if ((nlinhas[i] == 0) || (nlinhas[i] == 1))
			kLetter++;

		string letterK = takeLetter(kLetter);
		for (size_t j = 0; j < nlinhas[i]; j++)
			letterK += "\'";

		setsLetter[i] = letterK;
	}

}



void ChangeNames::generateOrderingGroupPoint(
	string fileName,
	std::vector<int> &uniqRcw,
	std::vector<int> &uniqCount,
	std::vector<int> &allNumbers,
	std::vector<std::string> &uniqPgroup)
{
	ifstream file_(fileName.c_str());
	allNumbers.resize(2);
	allNumbers[0] = 0;
	allNumbers[1] = 0;
	string firstLine;
	getline(file_, firstLine);
	string line;
	ReadWriteFormats rwf_;
	vector<int> allRcw;
	vector<string> allVgroup;
	vector<string> allPgroup;
	while (!file_.eof())
	{
		getline(file_, line);
		if (line == "")
			continue;
		stringstream convert;
		convert << line;
		int rcw;
		string vGroup, pGroup;
		int chiral, achiral;
		rwf_.takeRcwVgroupPointGroup(line, rcw, vGroup, pGroup);
		rwf_.takeRcwVgroupPointGroupNca(line, rcw, chiral, achiral, vGroup, pGroup);
		allNumbers[0] += chiral;
		allNumbers[1] += achiral;
		allRcw.push_back(rcw);
		allVgroup.push_back(vGroup);
		allPgroup.push_back(pGroup);
	}
	file_.close();

	AuxMath auxMath_;
	vector<int> instructions = auxMath_.vector_ordering(allRcw);
	auxMath_.vector_ordering_with_instructions(allVgroup, instructions);
	auxMath_.vector_ordering_with_instructions(allPgroup, instructions);

	// unique pGroup
	uniqRcw.push_back(allRcw[0]);
	uniqCount.push_back(1);
	uniqPgroup.push_back(allPgroup[0]);

	for (size_t i = 1; i < allRcw.size(); i++)
	{
		bool found = false;
		for (size_t j = 0; j < uniqPgroup.size(); j++)
		{
			if (allPgroup[i] == uniqPgroup[j])
			{
				if (allRcw[i] != uniqRcw[j])
				{
					cout << "point group problem: " << i << endl;
					exit(1);
				}
				uniqCount[j]++;
				found = true;
				break;
			}
		}
		if (!found)
		{
			uniqCount.push_back(1);
			uniqRcw.push_back(allRcw[i]);
			uniqPgroup.push_back(allPgroup[i]);
		}
	}

	//rwf_.symmetryGroupOrdering(uniqRcw, uniqPgroup, uniqCount);

	/*
	ofstream newCounting_("counting.csv");
	for (size_t i = 0; i < uniqpGroup.size(); i++)
	{
		newCounting_ << uniqpGroup[i] << " ;  ; ";
	}
	newCounting_ << endl;
	for (size_t i = 0; i < uniqCount.size(); i++)
	{
		newCounting_ << uniqCount[i] << " ; " << uniqRcw[i] << " ; ";
	}
	*/

}

void ChangeNames::generateOrderingGroupPoint(
	string fileName,
	std::vector<int> &uniqRcw,
	std::vector<int> &uniqCount,
	std::vector<int> &allNumbers,
	std::vector<std::string> &uniqPgroup,
	string &firstLine,
	std::vector<std::string> & listOfPermutations,
	std::vector<std::string> & listOfChiralities,
	std::vector<int> & listOfRcw,
	std::vector<std::string> & listOfPGroup)
{
	ifstream file_(fileName.c_str());
	allNumbers.resize(2);
	allNumbers[0] = 0;
	allNumbers[1] = 0;
	getline(file_, firstLine);
	string line;
	ReadWriteFormats rwf_;
	vector<int> allRcw;
	vector<string> allVgroup;
	vector<string> allPgroup;
	while (!file_.eof())
	{
		getline(file_, line);
		if (line == "")
			continue;
		stringstream convert;
		convert << line;
		int rcw;
		string vGroup, pGroup;
		int chiral, achiral;
		string permut;
		rwf_.takeAllElementsFromCode(line, rcw, chiral, achiral, vGroup, pGroup,permut);
		allNumbers[0] += chiral;
		allNumbers[1] += achiral;
		allRcw.push_back(rcw);
		allVgroup.push_back(vGroup);
		allPgroup.push_back(pGroup);
		listOfPermutations.push_back(permut);
		listOfRcw.push_back(rcw);
		listOfPGroup.push_back(pGroup);
		if (chiral == 1)
			listOfChiralities.push_back("c");
		else
			listOfChiralities.push_back("a");
	}
	file_.close();

	AuxMath auxMath_;
	vector<int> instructions = auxMath_.vector_ordering(allRcw);
	auxMath_.vector_ordering_with_instructions(allVgroup, instructions);
	auxMath_.vector_ordering_with_instructions(allPgroup, instructions);

	// unique pGroup
	uniqRcw.push_back(allRcw[0]);
	uniqCount.push_back(1);
	uniqPgroup.push_back(allPgroup[0]);

	for (size_t i = 1; i < allRcw.size(); i++)
	{
		bool found = false;
		for (size_t j = 0; j < uniqPgroup.size(); j++)
		{
			if (allPgroup[i] == uniqPgroup[j])
			{
				if (allRcw[i] != uniqRcw[j])
				{
					cout << "point group problem: " << i << endl;
					exit(1);
				}
				uniqCount[j]++;
				found = true;
				break;
			}
		}
		if (!found)
		{
			uniqCount.push_back(1);
			uniqRcw.push_back(allRcw[i]);
			uniqPgroup.push_back(allPgroup[i]);
		}
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

	vector<vultorGroup> auxVGroup = setVultorGroup(
		vultorGroupProbs,
		vultorGroupCounting,
		vultorGroupChirality);
	vGroup = auxVGroup;

	vector<double> entropyOrdering(auxVGroup.size());
	for (size_t i = 0; i < auxVGroup.size(); i++)
	{
		entropyOrdering[i] = auxVGroup[i].achiralN * auxVGroup[i].achiralProb
			+ auxVGroup[i].chiralN * auxVGroup[i].chiralProb;
	}
	AuxMath auxMath_;
	vector<int> instructions = auxMath_.vector_ordering(entropyOrdering);
	// ORDER VECTOR WITH INSTRUCTIONS
	for (size_t i = 0; i < instructions.size(); i += 2)
	{
		cout << "inversion" << endl;
		vGroup[instructions[i]] = auxVGroup[instructions[i + 1]];
		vGroup[instructions[i + 1]] = auxVGroup[instructions[i]];
	}
	for (size_t i = 0; i < vGroup.size(); i++)
	{
		vGroup[i].blockName[0] = takeLetter(i)[0];
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
		<< 2 * totalChiral + totalAchiral << " ; "
		<< totalChiral << " ; "
		<< totalAchiral << " ; ";


	int lastRce = vGroup[vGroup.size()-1].achiralN * vGroup[vGroup.size()-1].achiralProb
		+ vGroup[vGroup.size()-1].chiralN * vGroup[vGroup.size()-1].chiralProb;
	for (size_t i = 0; i < vGroup.size() - 1; i++)
	{
		int rce = vGroup[i].achiralN * vGroup[i].achiralProb + vGroup[i].chiralN * vGroup[i].chiralProb;
		if(rce % lastRce != 0)
			counting_ << fixed << setprecision(1) << (double)rce / (double)lastRce << ":";		
		else
			counting_ << rce / lastRce << ":";
	}
	counting_ << "1" << " ; ";

	for (size_t i = 0; i < vGroup.size(); i++)
	{
		int vGroupProb;
		if (vGroup[i].achiralProb == 0)
			vGroupProb = vGroup[i].chiralProb;
		else
			vGroupProb = vGroup[i].achiralProb;
		//int totalRce = 
		//	vGroup[i].chiralN * vGroup[i].chiralProb + 
		//	vGroup[i].achiralN * vGroup[i].achiralProb;

		counting_ << vGroup[i].chiralN << " ; "
			<< vGroup[i].achiralN << " ; "
			<< vGroupProb << " ; ";
			//<< totalRce << " ; ";

		/*
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
		*/
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


string ChangeNames::generateNewTypeLine(string pathRead, string &combination, int systemSize)
{
	string typeLineFileName = pathRead + "final-" + combination + "-atomTypes";
	ifstream typesFile_(typeLineFileName.c_str());
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
	int readingBid = -1;
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
		for (size_t i = 0; i < bidentates.size(); i++)
		{
			convert2 << bidentates[i] + 1 << "  ";
		}
	}
	return convert2.str();
}
