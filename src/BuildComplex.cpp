#include "BuildComplex.h"

#include <vector>
#include <iostream>

#include "AdjustSaParameters.h"
#include "ReadInput.h"
#include "Coordstructs.h"
#include "Ligand.h"
#include "MyExceptions.h"
#include "ComplexCreator.h"
#include "ControlMopac.h"

using namespace std;

BuildComplex::BuildComplex(){}

BuildComplex::~BuildComplex(){}

void BuildComplex::build()
{
	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.47976429810660104);

		ReadInput readInp_;
		if (!ReadLumpacViewInput(readInp_))
			return;

#ifdef _DEBUG
		readInp_.rePrintInput();
#endif

		bool terminate = buildLigands(readInp_);
		if (terminate) return;

		vector<CoordXYZ> allAtoms;
		if (!constructComplex(readInp_, saParameters_, allAtoms))
			return;

		runMopac(readInp_, allAtoms);
		return;
}

void BuildComplex::fitSA()
{
	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.47976429810660104);

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
		return sucess;
	}

	allAtoms = cpCreator.simulatedAnnealing();

	return sucess;
}

void BuildComplex::runMopac(ReadInput & readInp_, vector<CoordXYZ>& allAtoms)
{
	string mopacHeader = "RM1 BFGS PRECISE NOINTER XYZ T=10D GNORM=0.25 + \n NOLOG GEO-OK SCFCRT=1.D-10";
	string mopacFreq = "RM1 PRECISE NOINTER XYZ T=10D AUX THERMO FORCE + \n NOLOG GEO-OK SCFCRT=1.D-10";
	string mopacExecPath = "M2009_Ln_Orbitals.exe  ";
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



