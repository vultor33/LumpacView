#include "BuildComplex.h"

#include <vector>
#include <iostream>
#include <string>
#include <sstream>

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

void BuildComplex::build()
{
		ReadInput readInp_;
		if (!ReadLumpacViewInput(readInp_))
			return;

		build(readInp_);
}

void BuildComplex::build(string ligandName, int coordination, int charge, vector<string> options)
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
	readInp_.setProperties(options);

	build(readInp_);
}

void BuildComplex::build(ReadInput & readInp_)
{

	bool terminate = buildLigands(readInp_);
	if (terminate) return;

	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.5000000000000000);

	vector<CoordXYZ> allAtoms;
	if (!constructComplex(readInp_, saParameters_, allAtoms))
		return;

	runMopac2(readInp_, allAtoms);

}

void BuildComplex::checkIfIsSameIsomer(string xRayName)
{

	// PASSO 2 COMPARAR OS COMPLEXOS MONTADOS NOS DOIS CASOS
	RootMeanSquareDeviation rmsd_;
	vector<CoordXYZ> molXRay = rmsd_.readCoord(xRayName);
	for (size_t i = 0; i < molXRay.size(); i++)
	{



	}



















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
		double stretchDistance = 2.5e0;
		int saMaxIterations = 3000;
		ComplexCreator cpCreator(
			readInp_.allLigands,
			maxChelation,
			stretchDistance,
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
	double stretchDistance = 2.5e0;
	int saMaxIterations = 1000;
	ComplexCreator cpCreator(
		readInp_.allLigands,
		maxChelation,
		stretchDistance,
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
	string mopacExecPath = "M2009_Ln_Orbitals.exe";
	string mopacHeader = "RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	string mopacFreq = "RM1 PRECISE NOINTER XYZ T=10D AUX THERMO FORCE + \n NOLOG GEO-OK SCFCRT=1.D-10";
	ControlMopac controlMop(
		readInp_.getProjectName(),
		readInp_.getMetalName(),
		mopacHeader,
		mopacFreq,
		readInp_.getMetalParams(),
		mopacExecPath);
	bool sucess = controlMop.optimize(allAtoms);
	if (sucess) cout << "a estrutura otimizou com sucesso - tenha um bom dia" << endl;
	else {
		cout << "ocorreu algum problema na otimizacao do mopac, tente outra vez" << endl;
	}
}

void BuildComplex::runMopac2(ReadInput & readInp_, vector<CoordXYZ>& allAtoms)
{
	vector<string> options = readInp_.getOptions();

	string mopacExecPath = "M2009_Ln_Orbitals.exe";

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

bool BuildComplex::optimize(
	string mopacExecPath,
	vector<CoordXYZ> & allAtoms,
	vector<string> & options)
{
	vector<MopacParams> dummyParams;

	return optimize(mopacExecPath, allAtoms, options,dummyParams);
}


bool BuildComplex::optimize(
	string mopacExecPath,
	vector<CoordXYZ> & allAtoms,
	vector<string> & options,
	vector<MopacParams> & params)
{

	ReadQuantumOutput readmop_(options[0]);

	WriteQuantumInput writeMop_(options);

	string inputName = writeMop_.createInput(allAtoms);

	createParamsFile(options[3], params);

	system((mopacExecPath + " " + inputName).c_str());

	ReadQuantumOutput readGamess_(options[0]);

	readmop_.readOutput(inputName);

	allAtoms = readmop_.getCoordinates();

	if (allAtoms.size() == 0) return false;

	return true;
}

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


