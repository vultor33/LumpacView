#include "IsomersToMol.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "Coordstructs.h"
#include "Geometries.h"
#include "ReadWriteFormats.h"
#include "ChangeNames.h"

using namespace std;

IsomersToMol::IsomersToMol(){}

IsomersToMol::~IsomersToMol(){}

std::vector<std::string> IsomersToMol::readAllPermutations(
	std::string fileName,
	std::string fileFolder,
	vector<int> &atomTypes, 
	vector<int> &bidentateChosen)
{

	// read input
	ReadWriteFormats rwf_;
	size_t found = fileName.find("-");
	size_t foundSecond = fileName.find("-", found + 1, 1);
	size_t pointChar = fileName.find(".");
	string composition = fileName.substr(foundSecond + 1, pointChar - foundSecond - 1);
	int nBidentates;
	int coordination = rwf_.compositionToNumbers(composition, nBidentates);
	setParameters(coordination);
	ifstream fileIsomers_((fileFolder + fileName).c_str());

	rwf_.readAtomTypesAndBidentateChosenFile(
		fileIsomers_,
		atomTypes,
		bidentateChosen,
		coordination,
		nBidentates);
	string line;
	vector<string> allPermut;
	while (!fileIsomers_.eof())
	{
		vector<int> permutation = readCauchyNotationsEnantiomers(fileIsomers_);
		if (permutation.size() == 0)
			continue;
		string permtString = permutationToString0Correction(permutation);
		allPermut.push_back(permtString);
	}
	fileIsomers_.close();
	return allPermut;
}


//salvar codigo e grupo vultor dois dois
vector<int> IsomersToMol::findEnantiomerPair(
	std::string fileName,
	std::string fileFolder,
	vector<int> guessPermutation,
	vector<string> &pairCodes)
{
	for (size_t i = 0; i < guessPermutation.size(); i++)
		guessPermutation[i]++;
	vector<int> pair;
	ifstream fileIsomers_((fileFolder + fileName).c_str());
	string line;
	getline(fileIsomers_, line);
	vector<string> allPermut;
	bool first = true;
	vector<string> auxCode;
	vector<int> dummy;
	while (!fileIsomers_.eof())
	{
		vector<string> code1;
		vector<int> permutation = readCauchyNotationsEnantiomersAndTakeCode(fileIsomers_,code1);
		if ((first) && (permutation == guessPermutation))
		{
			vector<string> code2;
			pair = readCauchyNotationsEnantiomersAndTakeCode(fileIsomers_, code2);
			pairCodes.push_back(code1[0]);
			pairCodes.push_back(code1[1]);
			if (pair.size() != 0)
			{
				pairCodes.push_back(code2[0]);
				pairCodes.push_back(code2[1]);
			}
			break;
		}
		else if (first)
		{
			pair = permutation;
			auxCode = code1;
			first = false;
		}
		else if ((!first) && (permutation == guessPermutation))
		{
			pairCodes.push_back(code1[0]);
			pairCodes.push_back(code1[1]);
			pairCodes.push_back(auxCode[0]);
			pairCodes.push_back(auxCode[1]);
			break;
		}
		else if (permutation.size() == 0)
		{
			pair = dummy;
			first = true;
			continue;
		}
	}
	fileIsomers_.close();
	if (pair.size() != 0)
	{
		for (size_t i = 0; i < pair.size(); i++)
			pair[i]--;
	}
	return pair;
}

void IsomersToMol::printAllMolFromSpecifiedGeometry(
	int geoCode,
	std::string pathRead,
	std::string responseName)
{
	cout << "WARNING --- WORKS ONLY ON LINUX" << endl;

	Geometries geo_;
	ReadWriteFormats rwf_;
	string geomName = geo_.sizeToGeometryCode(geoCode);
	ifstream response_((pathRead + responseName).c_str());
	string line;
	string folder = geo_.sizeToGeometryCode(geoCode);

	system(("mkdir " + folder).c_str());
	while (!response_.eof())
	{
		getline(response_, line);
		cout << "line: " << line << endl;
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
		IsomersToMol ismol_;

		ismol_.printAllMol(allIsomersCombinationFile,
			pathRead,
			geoCode);
	}
}

// need geo code folder
void IsomersToMol::printAllMol(
	string fileName,
	string filePath,
	int geoCode)
{
	// read input
	size_t found = fileName.find("-");
	size_t foundSecond = fileName.find("-", found + 1, 1);
	size_t pointChar = fileName.find(".");
	string composition = fileName.substr(foundSecond + 1, pointChar - foundSecond - 1);
	int nBidentates;
	ReadWriteFormats rwf_;
	int coordination = rwf_.compositionToNumbers(composition, nBidentates);
	setParameters(coordination, geoCode);
	ifstream fileIsomers_((filePath + fileName).c_str());
	vector<int> atomTypes;
	vector<int> bidentateChosen;

	/* Labels - O, P, N, / bidentates  ---  format */
	rwf_.readAtomTypesAndBidentateChosenFileWithLabels(
		fileIsomers_,
		atomTypes,
		bidentateChosen,
		coordination,
		nBidentates);

	Geometries geo_;
	string folder = geo_.sizeToGeometryCode(geoCode);
	string folderCompositionZero = folder + "//" + composition;

	string folderComposition;

	string folderComp = composition;

	rwf_.ReplaceAll(folderComp,'(',"\\");
	rwf_.ReplaceAll(folderComp,')',"\\");
	//replace(folderComp.begin(),folderComp.end(),'(','[');
	//replace(folderComp.begin(),folderComp.end(),')',']');

	folderComposition = folder + "//" + folderComp;

	system(("mkdir " + folderComposition).c_str());

	rwf_.ReplaceAll(fileName,'(',"\\");
	rwf_.ReplaceAll(fileName,')',"\\");
	//replace(fileName.begin(),fileName.end(),'(','[');
	//replace(fileName.begin(),fileName.end(),')',']');

	system(("cp " + filePath + fileName + "  " + folderComposition).c_str());

	while (!fileIsomers_.eof())
	{
		string line;

		getline(fileIsomers_, line);
		if (line == "")
			continue;
		
		string vGroup, pGroup, vCode, setGroup, chirality;
		int rcw;
		vector<int> permutation;
		rwf_.takeAllElementsFromCodeNewSym(line, 
			coordination,
			rcw,
			chirality, 
			vGroup,
			pGroup,
			setGroup,
			permutation);

		ChangeNames chNames_;
		stringstream convertSymN;
		convertSymN << chNames_.symmetryNumberFromPointGroup(pGroup);
		string subFolderName = pGroup + "_" + chirality + "_" + convertSymN.str() + "_" + setGroup[0];

		system(("mkdir " + folderComposition + "//" + subFolderName).c_str());
		ofstream fileXyz_((folderCompositionZero + "//" + subFolderName + "//" + vGroup + ".mol2").c_str());

		for(size_t i = 0; i < permutation.size(); i++)
			permutation[i] = permutation[i] - 1;

		printMoleculeMolFormat(permutation, atomTypes, bidentateChosen, fileXyz_);

		fileXyz_.close();

	}
	fileIsomers_.close();
}

void IsomersToMol::printSingleMol(
	vector<int> &permutation, 
	vector<int> &atomTypes,
	vector<int> &bidentateChosen,
	string name)
{
	name += "-";
	vector<int> permutation2 = permutation;
	for (size_t i = 0; i < permutation2.size(); i++)
		permutation2[i] += 1;
	name += permutationToString(permutation2);	
	name += ".mol2";
	ofstream fileXyz_(name.c_str());
	printMoleculeMolFormat(permutation, atomTypes, bidentateChosen, fileXyz_);
	fileXyz_.close();

}


/* 
void IsomersToMol::readAtomTypesAndBidentateChosenFile(
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
		atomTypes.push_back(aux-1);
	}
	if (nBidentates != 0)
	{
		string dummy;
		auxline >> dummy;
		for(int i = 0; i < 2 * nBidentates; i++)
		{
			int aux;
			auxline >> aux;
			bidentateChosen.push_back(aux-1);
		}
	}
}
*/

vector<int> IsomersToMol::readCauchyNotationsEnantiomers(ifstream & openendFile_)
{
	int size = complex.size();
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
	}
	return notation;
}

vector<int> IsomersToMol::readCauchyNotationsEnantiomersAndTakeCode(
	ifstream & openendFile_,
	vector<string> &permutCodes)
{
	int size = complex.size();
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
	string vultorGroup = auxline.substr(brack2 + 2, 2);
	size_t key1 = auxline.find("{");
	string allCode = auxline.substr(key1, allCode.size() - key1);
	permutCodes.push_back(vultorGroup);
	permutCodes.push_back(allCode);

	stringstream line;
	line << permString;
	for (size_t i = 0; i < size; i++)
	{
		line >> notation[i];
	}
	return notation;
}



string IsomersToMol::permutationToString0Correction(vector<int> &permutation)
{
	stringstream permt;
	for (size_t i = 0; i < permutation.size(); i++)
		permt << (permutation[i] - 1) << " ";
	return permt.str();
}

string IsomersToMol::permutationToString(vector<int> &permutation)
{
	stringstream permt;
	for (size_t i = 0; i < permutation.size(); i++)
		permt << permutation[i] << " ";
	return permt.str();
}

void IsomersToMol::printMoleculeMolFormat(
	vector<int> & permutation,
	vector<int> & atomTypes,
	const vector<int> & bidentateAtomsChosen,
	ofstream & printFile_)
{

	double scale = 3.0e0;
	vector<string> atoms(permutation.size());
	for (size_t i = 0; i < atomTypes.size(); i++)
		atoms[i] = atomLabels[atomTypes[i]];

	if (permutation.size() == 0)
	{
		permutation.resize(atoms.size());
		for (size_t i = 0; i < atoms.size(); i++)
			permutation[i] = i;
	}
	//defining bidentate bonds
	vector<int> bidentateAtomsChosenRotated = applyPermutationCoordinates(permutation, atoms, bidentateAtomsChosen);
	int nAtoms = complex.size();
	int nBonds = nAtoms + (bidentateAtomsChosenRotated.size() / 2);
	printFile_ << "#   generated by: ComplexBuild" << endl << endl << endl;
	printFile_ << "@<TRIPOS>MOLECULE" << endl;
	printFile_ << "name" << endl;
	printFile_ << nAtoms + 1 << "  " << nBonds << endl;
	printFile_ << "SMALL " << endl << "NO_CHARGES" << endl << endl << endl;
	printFile_ << "@<TRIPOS>ATOM" << endl;

	//lantanide
	double zero = 0.0e0;
	printFile_ << " 1 A1  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		//<< fixed << setprecision(4) << setw(8) << 0.0002e0 << "  "
		//<< fixed << setprecision(4) << setw(8) << 0.0487e0 << "  "
		<< "Au " << endl;


	for (int i = 0; i < nAtoms; i++)
	{
		printFile_ << " " << i + 2 << " A" << i + 2 << "  "
			<< fixed << setprecision(4) << setw(8) << complex[i].x * scale << "  "
			<< fixed << setprecision(4) << setw(8) << complex[i].y * scale << "  "
			<< fixed << setprecision(4) << setw(8) << complex[i].z * scale << "  "
			<< complex[i].atomlabel << endl;
	}
	printFile_ << "@<TRIPOS>BOND" << endl;
	int k = 1;
	for (size_t i = 0; i < nAtoms; i++)
	{
		printFile_ << " " << k << " 1  " << i + 2 << " 1 " << endl;
		k++;
	}
	for (size_t i = 0; i < bidentateAtomsChosen.size(); i += 2)
	{
		printFile_ << " " << k << " " << bidentateAtomsChosenRotated[i] + 2 << "  "
			<< bidentateAtomsChosenRotated[i + 1] + 2
			<< " 1 " << endl;
		k++;
	}
	printFile_ << "M  END" << endl;
}

std::vector<int> IsomersToMol::applyPermutationCoordinates(
	const std::vector<int> & permutation,
	const std::vector<std::string> & atoms,
	const std::vector<int> & bidentateAtomsChosen)
{
	size_t size = atoms.size();
	for (size_t i = 0; i < size; i++)
	{
		complex[i].atomlabel = atoms[permutation[i]];
	}
	vector<int> bidentateAtomsChosenRotated = bidentateAtomsChosen;
	for (size_t i = 0; i < bidentateAtomsChosenRotated.size(); i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			if (permutation[j] == bidentateAtomsChosenRotated[i])
			{
				bidentateAtomsChosenRotated[i] = j;
				break;
			}
		}
	}
	return bidentateAtomsChosenRotated;
}

/* fredapagar
int IsomersToMol::stringToNumber(string entryString, int & nBidentates)
{
	// cut ligands codes
	vector<string> allCodes;
	int coordination = 0;
	nBidentates = 0;
	bool bidentate;
	for (size_t i = 1; i < entryString.size();i++)
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
		if ((i == entryString.size())||
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

	coordination = 0;
	nBidentates = 0;
	for (size_t i = 0; i < typeCode1.size(); i++)
		coordination += typeCode1[i];
	for (size_t i = 0; i < typeCode2.size(); i++)
	{
		coordination += 2 * typeCode2[i];
		nBidentates += typeCode2[i];
	}
	for (size_t i = 0; i < typeCode3.size(); i++)
	{
		coordination += 2 * typeCode3[i];
		nBidentates += typeCode3[i];
	}

	return coordination;
}

int IsomersToMol::codeToType(string code)
{
	if (code[0] == 'M')
		return 0;
	else if (code[1] == 's')
		return 1;
	else if (code[1] == 'a')
		return 2;

	return -1;
}

void IsomersToMol::addEqual(
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

void IsomersToMol::addDifferent(
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
*/


string IsomersToMol::getAtomLabelI(int I)
{
	return atomLabels[I];
}

// molekel parameters
void IsomersToMol::setParameters(
	int coordination, 
	int geoCode)
{
	//molekel
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

	Geometries geo_;
	double cutAngle;
	vector<int> reflection;
	geo_.selectGeometry(geoCode, complex, cutAngle, reflection);

}

// raswin parameters
void IsomersToMol::setParameters(
	int coordination)
{
	//cout << "IsomersToMol::setParameters need attention - RASWIN format" << endl;
	//exit(1);

	atomLabels.resize(12);
	atomLabels[0] = "Ca";
	atomLabels[1] = "O";
	atomLabels[2] = "H";
	atomLabels[3] = "Mg";
	atomLabels[4] = "Be";
	atomLabels[5] = "I";
	atomLabels[6] = "Cl";
	atomLabels[7] = "Br";
	atomLabels[8] = "Na";
	atomLabels[9] = "He";
	atomLabels[10] = "N";
	atomLabels[11] = "C";

	Geometries geo_;
	double cutAngle;
	vector<int> reflection;

	switch (coordination)
	{
	case 4:
		geo_.selectGeometry(40, complex, cutAngle, reflection);
		break;

	case 5:
		geo_.selectGeometry(51, complex, cutAngle, reflection);
		break;

	case 6:
		geo_.selectGeometry(60, complex, cutAngle, reflection);
		break;

	case 7:
		geo_.selectGeometry(70, complex, cutAngle, reflection);
		break;

	case 8:
		geo_.selectGeometry(81, complex, cutAngle, reflection);
		break;

	case 9:
		geo_.selectGeometry(90, complex, cutAngle, reflection);
		break;

	case 10:
		geo_.selectGeometry(100, complex, cutAngle, reflection);
		break;

	case 11:
		geo_.selectGeometry(110, complex, cutAngle, reflection);
		break;

	case 12:
		geo_.selectGeometry(120, complex, cutAngle, reflection);
		break;

	default:
		cout << "Coordination not found, check input for errors" << endl;
		exit(1);
		break;
	}
}


void IsomersToMol::printCoordXyz(vector<CoordXYZ> &coord)
{
	ofstream xyz_("teste.xyz");
	xyz_ << coord.size() << endl
		<< "titulo" << endl;
	for (size_t i = 0; i < coord.size(); i++)
	{
		xyz_ << "H " << "  "
			<< coord[i].x << "  "
			<< coord[i].y << "  "
			<< coord[i].z << endl;
	}


}
