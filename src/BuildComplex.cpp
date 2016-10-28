#include "BuildComplex.h"

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

#include "AdjustSaParameters.h"
#include "ReadInput.h"
#include "Coordstructs.h"
#include "Ligand.h"
#include "MyExceptions.h"
#include "ComplexCreator.h"
#include "ControlMopac.h"
#include "RootMeanSquareDeviation.h"
#include "WriteQuantumInput.h"
#include "ReadQuantumOutput.h"

using namespace std;

BuildComplex::BuildComplex(){}

BuildComplex::~BuildComplex(){}

vector<CoordXYZ> BuildComplex::build()
{
	ReadInput readInp_ = activateReadInput();
	return build(readInp_);
}

vector<CoordXYZ> BuildComplex::buildCompletingWithWater(
	string ligandName, 
	int coordination, 
	int charge, 
	vector<string> options, 
	string mopacExecPath,
	int & nLigandAtoms)
{
	Ligand mol;
	mol.setLigandCoordinates(ligandName);

	int nCl = 3 - charge;

	int nH2o = coordination - mol.getChelation() - nCl;

	vector<Ligand> allLigands(1 + nCl + nH2o);

	Ligand Cl_, H2o_;
	
	setClandH2oToCompleteCoordination(Cl_, H2o_);

	allLigands[0] = mol;
	int k = 1;
	for (int i = 0; i < nCl; i++)
	{
		allLigands[k] = Cl_;
		k++;
	}
	for (int i = 0; i < nH2o; i++)
	{
		allLigands[k] = H2o_;
		k++;
	}

	ReadInput readInp_;
	readInp_.allLigands = allLigands;
	readInp_.setProperties(options, mopacExecPath);

	nLigandAtoms = mol.getNatoms();

	return build(readInp_);
}

vector<CoordXYZ> BuildComplex::build(ReadInput & readInp_)
{
	bool terminate = buildLigands(readInp_);
	if (terminate) return vector<CoordXYZ>();

	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.5000000000000000);

	vector<CoordXYZ> allAtoms;
	if (!constructComplex(readInp_, saParameters_, allAtoms))
		return allAtoms;

	runMopac(readInp_, allAtoms); //allAtoms modified - Metal always on origin

	return allAtoms;
}

void BuildComplex::makeComplexOptimizingInMopac(string ligandName, int coordination, int charge, vector<string> options, string mopacExecPath)
{
	int nAtoms;
	vector<CoordXYZ> allAtoms = buildCompletingWithWater(ligandName, coordination, charge, options, mopacExecPath, nAtoms);
	vector<CoordXYZ> ligandCreated(nAtoms);
	for (int i = 0; i < nAtoms; i++)
		ligandCreated[i] = allAtoms[i + 1];

	Ligand newLigand_;
	newLigand_.setLigandCoordinates(ligandName);
	newLigand_.setNewCoordinates(ligandCreated);
	double metalLigandDistance = newLigand_.distanceX1ToPoint(allAtoms[0].x, allAtoms[0].y, allAtoms[0].z);
	newLigand_.initializeLigand();

	string newLigandName = "Lumpac-View-Ligand-" + ligandName;

	ofstream newLigFile_(newLigandName.c_str());
	newLigFile_ << newLigand_.getNatoms() << endl;
	int chelation = newLigand_.getChelation();
	string chelationName;
	if (chelation == 1)
		chelationName = "monodentate";
	else if (chelation == 2)
		chelationName = "bidentate";
	else
		chelationName = "tridentate";
	
	newLigFile_ << newLigand_.getFormalName() << " "
		<< chelationName << " "
		<< metalLigandDistance << endl;

	newLigand_.printLigand(newLigFile_);
	newLigFile_.close();
}

vector<Ligand> BuildComplex::assembleComplexWithoutSA()
{

	vector<int> ligandsPermutation; 
	vector<string>  inputInformations;
	return assembleComplexWithoutSA(ligandsPermutation,inputInformations);
}

vector<Ligand> BuildComplex::assembleComplexWithoutSA(vector<int> & ligandsPermutation)
{
	vector<string>  inputInformations;
	return assembleComplexWithoutSA(ligandsPermutation,inputInformations);
}

vector<Ligand> BuildComplex::assembleComplexWithoutSA(vector<int> & ligandsPermutation, vector<string> & inputInformations)
{
	ReadInput readInp_ = activateReadInput(inputInformations);

	bool terminate = buildLigands(readInp_);

	if (ligandsPermutation.size() == 0)
	{
		int allChelation = 0;
		for (size_t i = 0; i < readInp_.allLigands.size(); i++)
		{
			allChelation += readInp_.allLigands[i].getChelation();
		}
		actualLigandPermutation.resize(allChelation);
		for (int i = 0; i < allChelation; i++)
		{
			actualLigandPermutation[i] = i;
		}
	}
	else
		actualLigandPermutation = ligandsPermutation;

	if (terminate) return vector<Ligand>();

	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.5000000000000000);

	vector<CoordXYZ> allAtoms;

	int maxChelation = 12;
	int saMaxIterations = 1000;
	ComplexCreator cpCreator(
		readInp_.allLigands,
		maxChelation,
		saMaxIterations,
		saParameters_.maxAlfaAngle,
		saParameters_.maxBetaAngle,
		saParameters_.saTemperatureUpdate,
		saParameters_.saInitialTemperature,
		saParameters_.saAcceptance);
	bool sucess = cpCreator.start(actualLigandPermutation);
	vector<Ligand> allLigands = cpCreator.getLigandsCreated();

	return allLigands;
}

int BuildComplex::getLigandsNumber()
{
	ReadInput readInp_ = activateReadInput();
	return readInp_.allLigands.size();
}


ReadInput BuildComplex::activateReadInput()
{
	vector<string> inputInformations;
	activateReadInput(inputInformations);
}

ReadInput BuildComplex::activateReadInput(vector< string > & inputInformations)
{
	ReadInput readInp_;
	if (inputInformations.size() == 0)
	{
		if (!ReadLumpacViewInput(readInp_))
		{
			cout << "error on input" << endl;
			exit(1);
		}
	}
	else
		readInp_.setLumpacViewInputWithoutFile(inputInformations);

	vector<string> options = readInp_.getOptions();
	options[0] = "mopac2009";
	options[2] = " CHARGE=1 NOLOG GEO-OK SCFCRT=1.D-10";
	readInp_.setProperties(options, readInp_.getMopacExecPath());
	readInp_.setProperties(readInp_.getOptions(), "M2009_Ln_Orbitals.exe");
	return readInp_;
}


void BuildComplex::fitSA()
{
	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.5000000000000000);

	bool paramLimits = saParameters_.takeParametersFromFile();
	if (!paramLimits)
	{
		remove("fitness.txt");
		ofstream fit_("fitness.txt");
		fit_ << 5000 << endl;
		fit_.close();
	}
	else
	{
		ReadInput readInp_;
		if (!ReadLumpacViewInput(readInp_))
			return;

#ifdef _DEBUG
		readInp_.rePrintInput();
#endif

		bool terminate = buildLigands(readInp_);
		if (terminate) return;

		vector<CoordXYZ> allAtoms;

		int maxChelation = 10;
		int saMaxIterations = 3000;
		ComplexCreator cpCreator(
			readInp_.allLigands,
			maxChelation,
			saMaxIterations,
			saParameters_.maxAlfaAngle,
			saParameters_.maxBetaAngle,
			saParameters_.saTemperatureUpdate,
			saParameters_.saInitialTemperature,
			saParameters_.saAcceptance);

		bool sucess = cpCreator.start();
		if (sucess) cout << "complexo iniciado com sucesso" << endl;
		else {
			cout << "um problema aconteceu nos ligantes" << endl;
			return;
		}

		allAtoms = cpCreator.simulatedAnnealing();

		remove("fitness.txt");
		ofstream fit_("fitness.txt");
		fit_ << cpCreator.getIterationOfLowestFit() << endl;
		fit_.close();
	}
	return;
}


bool BuildComplex::ReadLumpacViewInput(ReadInput & readInp_)
{
	try {
		readInp_.readLumpacViewInput();
	}
	catch (MyExceptions& e) {
		cout << e.what() << endl;
		return false;
	}
	catch (...) {
		cout << "unknown error on input - check it or contact developers" << endl;
		return false;
	}
	return true;
}

bool BuildComplex::buildLigands(ReadInput & readInp_)
{
	bool sucess;
	bool terminateExecution = false;
	for (size_t i = 0; i < readInp_.allLigands.size(); i++) {
		sucess = readInp_.allLigands[i].initializeLigand();
		if (!sucess) {
			cout << "Can't build ligand:  " << i << " check input" << endl;
			terminateExecution = true;
		}
	}
	return terminateExecution;
}

bool BuildComplex::constructComplex(ReadInput & readInp_, const AdjustSaParameters & saParameters_, vector<CoordXYZ>& allAtoms)
{
	int maxChelation = 12;
	int saMaxIterations = 1000;
	ComplexCreator cpCreator(
		readInp_.allLigands,
		maxChelation,
		saMaxIterations,
		saParameters_.maxAlfaAngle,
		saParameters_.maxBetaAngle,
		saParameters_.saTemperatureUpdate,
		saParameters_.saInitialTemperature,
		saParameters_.saAcceptance);

	bool sucess = cpCreator.start();
	if (sucess) cout << "complexo iniciado com sucesso" << endl;
	else {
		cout << "um problema aconteceu nos ligantes" << endl;
		return sucess;
	}

	allAtoms = cpCreator.simulatedAnnealing();

	return sucess;
}

void BuildComplex::runMopac(ReadInput & readInp_, vector<CoordXYZ>& allAtoms)
{
	vector<string> options = readInp_.getOptions();

	string mopacExecPath = readInp_.getMopacExecPath();

	optimize(mopacExecPath, allAtoms, options);

}

void BuildComplex::setClandH2oToCompleteCoordination(Ligand & Cl_, Ligand & H2o_)
{
	vector<CoordXYZ> cl(1);
	cl[0].atomlabel = "Cl";
	cl[0].x = 0.0e0;
	cl[0].y = 0.0e0;
	cl[0].z = 0.0e0;
	string titleInfoCl = "Cl monodentate 2";
	Cl_.setLigandCoordinates(cl, titleInfoCl);

	vector<CoordXYZ> h2o(3);
	h2o[0].atomlabel = "O";
	h2o[0].x = 0.0e0;
	h2o[0].y = 0.0e0;
	h2o[0].z = 0.0e0;
	h2o[1].atomlabel = "H";
	h2o[1].x = 0.757e0;
	h2o[1].y = 0.586e0;
	h2o[1].z = 0.0e0;
	h2o[2].atomlabel = "H";
	h2o[2].x = -0.757e0;
	h2o[2].y = 0.586e0;
	h2o[2].z = 0.0e0;
	string titleInfoh2o = "H2O monodentate 2";
	H2o_.setLigandCoordinates(h2o, titleInfoh2o);
}

void BuildComplex::runMopac(vector<string> options, string mopacExecPath, vector<CoordXYZ>& allAtoms)
{
	optimize(mopacExecPath, allAtoms, options);
};

void BuildComplex::runMopacAndPrint(vector<string> options, string mopacExecPath, vector<CoordXYZ>& allAtoms)
{
	runMopac(options, mopacExecPath, allAtoms);

	ReadQuantumOutput readmop_(options[0]);

	readmop_.readOutput(options[1]);

	allAtoms = readmop_.getCoordinates();

	string fName = options[1] + ".xyz";
	
	ofstream pr_;
	pr_.open(fName.c_str(), std::ofstream::out | std::ofstream::app);
	pr_ << allAtoms.size() << endl << "E:  " << setprecision(16) << readmop_.getEnergy() << endl;
	for (size_t i = 0; i < allAtoms.size(); i++)
	{
		pr_ << allAtoms[i].atomlabel << "  "
			<< allAtoms[i].x << "  "
			<< allAtoms[i].y << "  "
			<< allAtoms[i].z << endl;
	}
	pr_.close();
};


bool BuildComplex::optimize(
	string mopacExecPath,
	vector<CoordXYZ> & allAtoms,
	vector<string> & options,
	vector<MopacParams> params)
{
	WriteQuantumInput writeMop_(options);

	string inputName = writeMop_.createInput(allAtoms);

	createParamsFile(options[3], params);

	system((mopacExecPath + " " + inputName).c_str());

	return true;
};


void BuildComplex::createParamsFile(string paramsName, vector<MopacParams>& params)
{
	if (params.size() == 0)
		return;

	string name = paramsName + ".inp";
	remove(name.c_str());
	ofstream paramsFile_(name.c_str());

	for (size_t i = 0; i < params.size(); i++)
		paramsFile_ << params[i].paramName << "  " << params[i].paramValue << endl;

	paramsFile_ << "  END     " << endl;
	paramsFile_.close();
}


