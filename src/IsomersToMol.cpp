#include "IsomersToMol.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "Coordstructs.h"

using namespace std;

IsomersToMol::IsomersToMol(){}

IsomersToMol::~IsomersToMol(){}

std::vector<std::string> IsomersToMol::readAllPermutations(
	std::string fileName, 
	vector<int> &atomTypes, 
	vector<int> &bidentateChosen)
{
	// read input
	size_t found = fileName.find("-");
	size_t foundSecond = fileName.find("-", found + 1, 1);
	size_t pointChar = fileName.find(".");
	string composition = fileName.substr(foundSecond + 1, pointChar - foundSecond - 1);
	int nBidentates;
	int coordination = stringToNumber(composition, nBidentates);
	setParameters(coordination);
	ifstream fileIsomers_(fileName.c_str());

	readAtomTypesAndBidentateChosenFile(
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

vector<int> IsomersToMol::findEnantiomerPair(
	std::string fileName,
	vector<int> guessPermutation)
{
	for (size_t i = 0; i < guessPermutation.size(); i++)
		guessPermutation[i]++;
	vector<int> pair;
	ifstream fileIsomers_(fileName.c_str());
	string line;
	getline(fileIsomers_, line);
	vector<string> allPermut;
	bool first = true;
	while (!fileIsomers_.eof())
	{
		vector<int> permutation = readCauchyNotationsEnantiomers(fileIsomers_);
		if ((first) && (permutation == guessPermutation))
		{
			pair = readCauchyNotationsEnantiomers(fileIsomers_);
			break;
		}
		else if (first)
		{
			pair = permutation;
			first = false;
		}
		else if ((!first) && (permutation == guessPermutation))
		{
			break;
		}
		else if (permutation.size() == 0)
		{
			first = true;
			continue;
		}
	}
	fileIsomers_.close();
	for (size_t i = 0; i < pair.size(); i++)
		pair[i]--;
	return pair;
}


void IsomersToMol::printAllMol(string fileName)
{

	// read input
	size_t found = fileName.find("-");
	size_t foundSecond = fileName.find("-", found + 1, 1);
	size_t pointChar = fileName.find(".");
	string composition = fileName.substr(foundSecond + 1, pointChar - foundSecond - 1);
	int nBidentates;
	int coordination = stringToNumber(composition, nBidentates);
	setParameters(coordination);
	ifstream fileIsomers_(fileName.c_str());

	vector<int> atomTypes;
	vector<int> bidentateChosen;
	readAtomTypesAndBidentateChosenFile(
		fileIsomers_, 
		atomTypes, 
		bidentateChosen, 
		coordination,
		nBidentates);
	string line;
	while (!fileIsomers_.eof())
	{
		vector<int> permutation = readCauchyNotationsEnantiomers(fileIsomers_);
		if (permutation.size() == 0)
			continue;
		string permtString = permutationToString0Correction(permutation);
		ofstream fileXyz_((composition + "-" + permtString + ".mol2").c_str());
		printMoleculeMolFormat(permutation, atomTypes, bidentateChosen, fileXyz_);
		fileXyz_.close();
	}
	fileIsomers_.close();

	cout << "FileToIsomers ended normally" << endl
		<< "Press enter to exit" << endl;
	string dummy;
	getline(cin, dummy);
}

void IsomersToMol::printSingleMol(
	vector<int> &permutation, 
	vector<int> &atomTypes,
	vector<int> &bidentateChosen,
	string name)
{
	name += "-";
	name += permutationToString(permutation);
	name += ".mol2";
	ofstream fileXyz_(name.c_str());
	printMoleculeMolFormat(permutation, atomTypes, bidentateChosen, fileXyz_);
	fileXyz_.close();

}


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

	size_t brack1 = auxline.find("[");
	size_t brack2 = auxline.find("]");
	string permString = auxline.substr(brack1 + 1, brack2 - brack1 - 1);
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
	printFile_ << " 1 A1 "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< fixed << setprecision(4) << setw(8) << zero << "  "
		<< "Au " << endl;

	for (int i = 0; i < nAtoms; i++)
	{
		printFile_ << " " << i + 2 << " A" << i + 1 << "  "
			<< fixed << setprecision(4) << setw(8) << complex[i].x * scale << "  "
			<< fixed << setprecision(4) << setw(8) << complex[i].y * scale << "  "
			<< fixed << setprecision(4) << setw(8) << complex[i].z * scale << "  "
			<< complex[i].atomlabel << endl;
	}
	printFile_ << "@<TRIPOS>BOND" << endl;
	int k = 1;
	for (size_t i = 0; i < nAtoms; i++)
	{
		printFile_ << " " << k << " 1  " << i + 2 << endl;
		k++;
	}
	for (size_t i = 0; i < bidentateAtomsChosen.size(); i += 2)
	{
		printFile_ << " " << k << " " << bidentateAtomsChosenRotated[i] + 2 << "  "
			<< bidentateAtomsChosenRotated[i + 1] + 2
			<< endl;
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

int IsomersToMol::stringToNumber(string entryString, int & nBidentates)
{
	// cut ligands codes
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

	int coordination = 0;
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


void IsomersToMol::setParameters(int coordination)
{
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

	switch (coordination)
	{
	case 4:
		complex.resize(4);
		complex[0].x = 0.0000000;
		complex[0].y = 0.0000000;
		complex[0].z = 1.0000000;
		complex[1].x = 0.81649658;
		complex[1].y = 0.47140452;
		complex[1].z = -0.33333333;
		complex[2].x = -0.81649658;
		complex[2].y = 0.47140452;
		complex[2].z = -0.33333333;
		complex[3].x = 0.00000000;
		complex[3].y = -0.94280904;
		complex[3].z = -0.33333333;
		break;

	case 5:
		complex.resize(5);
		complex[0].x = 0.000e0;
		complex[0].y = 0.000e0;
		complex[0].z = 1.000e0;
		complex[1].x = 1.000e0;
		complex[1].y = 0.000e0;
		complex[1].z = 0.000e0;
		complex[2].x = -1.000e0;
		complex[2].y = 0.000e0;
		complex[2].z = 0.000e0;
		complex[3].x = 0.000e0;
		complex[3].y = 0.86602540;
		complex[3].z = -0.50000000;
		complex[4].x = 0.000e0;
		complex[4].y = -0.86602540;
		complex[4].z = -0.50000000;
		break;

	case 6:
		complex.resize(6);
		complex[0].x = 0.00000000;
		complex[0].y = 0.00000000;
		complex[0].z = 1.00000000;
		complex[1].x = 1.00000000;
		complex[1].y = 0.00000000;
		complex[1].z = 0.00000000;
		complex[2].x = 0.00000000;
		complex[2].y = 1.00000000;
		complex[2].z = 0.00000000;
		complex[3].x = -1.00000000;
		complex[3].y = 0.00000000;
		complex[3].z = 0.00000000;
		complex[4].x = 0.00000000;
		complex[4].y = -1.00000000;
		complex[4].z = 0.00000000;
		complex[5].x = 0.00000000;
		complex[5].y = 0.00000000;
		complex[5].z = -1.00000000;
		break;

	case 7:
		complex.resize(7);
		complex[0].x = 0.00000000;
		complex[0].y = 0.00000000;
		complex[0].z = 1.00000000;
		complex[1].x = 0.97767167;
		complex[1].y = 0.00000000;
		complex[1].z = 0.21013831;
		complex[2].x = 0.16977090;
		complex[2].y = 0.96281864;
		complex[2].z = 0.21013831;
		complex[3].x = -0.91871085;
		complex[3].y = 0.33438340;
		complex[3].z = 0.21013831;
		complex[4].x = -0.48883583;
		complex[4].y = -0.84668850;
		complex[4].z = 0.21013831;
		complex[5].x = 0.36282725;
		complex[5].y = -0.62843523;
		complex[5].z = -0.68805926;
		complex[6].x = -0.26010411;
		complex[6].y = 0.45051354;
		complex[6].z = -0.85403946;
		break;

	case 8:
		complex.resize(8);
		complex[0].x = 0.000000;
		complex[0].y = 0.000000;
		complex[0].z = 1.00000000;
		complex[1].x = 0.96528366;
		complex[1].y = 0.000000;
		complex[1].z = 0.26120388;
		complex[2].x = -0.56545007;
		complex[2].y = 0.78232905;
		complex[2].z = 0.26120388;
		complex[3].x = -0.88247541;
		complex[3].y = -0.39116453;
		complex[3].z = 0.26120387;
		complex[4].x = 0.19991679;
		complex[4].y = -0.94435471;
		complex[4].z = 0.26120387;
		complex[5].x = 0.39983358;
		complex[5].y = 0.78232905;
		complex[5].z = -0.47759225;
		complex[6].x = -0.59975037;
		complex[6].y = 0.16202565;
		complex[6].z = -0.78361162;
		complex[7].x = 0.48264183;
		complex[7].y = -0.39116453;
		complex[7].z = -0.78361162;
		break;

	case 9:
		complex.resize(9);
		complex[0].x = 0.000000;
		complex[0].y = 0.000000;
		complex[0].z = 1.000000;
		complex[1].x = -0.23570226;
		complex[1].y = 0.91287093;
		complex[1].z = 0.33333333;
		complex[2].x = -0.94280904;
		complex[2].y = 0.00000000;
		complex[2].z = 0.33333333;
		complex[3].x = 0.23570226;
		complex[3].y = -0.91287093;
		complex[3].z = 0.33333333;
		complex[4].x = 0.94280904;
		complex[4].y = 0.00000000;
		complex[4].z = 0.33333333;
		complex[5].x = 0.53033009;
		complex[5].y = 0.68465320;
		complex[5].z = -0.50000000;
		complex[6].x = -0.53033009;
		complex[6].y = -0.68465320;
		complex[6].z = -0.50000000;
		complex[7].x = -0.58925565;
		complex[7].y = 0.45643546;
		complex[7].z = -0.66666667;
		complex[8].x = 0.58925565;
		complex[8].y = -0.45643546;
		complex[8].z = -0.66666667;
		break;

	case 10:
		complex.resize(10);
		complex[0].x = 0.000000;
		complex[0].y = 0.000000;
		complex[0].z = 1.000000;
		complex[1].x = 0.91458473;
		complex[1].y = 0.00000000;
		complex[1].z = 0.40439433;
		complex[2].x = 0.26335401;
		complex[2].y = 0.87584810;
		complex[2].z = 0.40439432;
		complex[3].x = -0.76291954;
		complex[3].y = 0.50439965;
		complex[3].z = 0.40439433;
		complex[4].x = -0.38787671;
		complex[4].y = -0.82826136;
		complex[4].z = 0.40439433;
		complex[5].x = 0.52670802;
		complex[5].y = -0.82826136;
		complex[5].z = -0.19121135;
		complex[6].x = -0.89351054;
		complex[6].y = -0.25145533;
		complex[6].z = -0.37203377;
		complex[7].x = 0.67837321;
		complex[7].y = 0.50439965;
		complex[7].z = -0.53421979;
		complex[8].x = -0.39464230;
		complex[8].y = 0.69067213;
		complex[8].z = -0.60599460;
		complex[9].x = 0.02107419;
		complex[9].y = -0.25145533;
		complex[9].z = -0.96763944;
		break;

	case 11:
		complex.resize(11);
		complex[0].x = 0.000000;
		complex[0].y = 0.000000;
		complex[0].z = 1.000000;
		complex[1].x = 0.89442719;
		complex[1].y = 0.00000000;
		complex[1].z = 0.44721360;
		complex[2].x = 0.27639320;
		complex[2].y = 0.85065081;
		complex[2].z = 0.44721360;
		complex[3].x = -0.72360680;
		complex[3].y = 0.52573111;
		complex[3].z = 0.44721360;
		complex[4].x = -0.72360680;
		complex[4].y = -0.52573111;
		complex[4].z = 0.44721360;
		complex[5].x = 0.27639320;
		complex[5].y = -0.85065081;
		complex[5].z = 0.44721360;
		complex[6].x = 0.72360680;
		complex[6].y = 0.52573111;
		complex[6].z = -0.44721360;
		complex[7].x = -0.27639320;
		complex[7].y = 0.85065081;
		complex[7].z = -0.44721360;
		complex[8].x = -0.89442719;
		complex[8].y = 0.00000000;
		complex[8].z = -0.44721360;
		complex[9].x = -0.27639320;
		complex[9].y = -0.85065081;
		complex[9].z = -0.44721360;
		complex[10].x = 0.00000000;
		complex[10].y = 0.00000000;
		complex[10].z = -1.00000000;
		break;

	case 12:
		complex.resize(12);
		complex[0].x = 0.00000000;
		complex[0].y = 0.00000000;
		complex[0].z = 1.00000000;
		complex[1].x = 0.89442719;
		complex[1].y = 0.00000000;
		complex[1].z = 0.44721360;
		complex[2].x = 0.27639320;
		complex[2].y = 0.85065081;
		complex[2].z = 0.44721360;
		complex[3].x = -0.72360680;
		complex[3].y = 0.52573111;
		complex[3].z = 0.44721360;
		complex[4].x = -0.72360680;
		complex[4].y = -0.52573111;
		complex[4].z = 0.44721360;
		complex[5].x = 0.27639320;
		complex[5].y = -0.85065081;
		complex[5].z = 0.44721360;
		complex[6].x = 0.72360680;
		complex[6].y = 0.52573111;
		complex[6].z = -0.44721360;
		complex[7].x = -0.27639320;
		complex[7].y = 0.85065081;
		complex[7].z = -0.44721360;
		complex[8].x = -0.89442719;
		complex[8].y = 0.00000000;
		complex[8].z = -0.44721360;
		complex[9].x = -0.27639320;
		complex[9].y = -0.85065081;
		complex[9].z = -0.44721360;
		complex[10].x = 0.72360680;
		complex[10].y = -0.52573111;
		complex[10].z = -0.44721360;
		complex[11].x = 0.00000000;
		complex[11].y = 0.00000000;
		complex[11].z = -1.00000000;
		break;

	default:
		cout << "Coordination not found, check input for errors" << endl;
		exit(1);
		break;
	}
}
