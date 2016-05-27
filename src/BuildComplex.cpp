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
#include "RootMeanSquareDeviation.h"

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

	string projectName = options[0];
	string metalName = options[1];
	string metalParams = options[2];
	readInp_.setProperties(projectName, metalName, metalParams);
}

void BuildComplex::build(ReadInput & readInp_)
{

#ifdef _DEBUG
	readInp_.rePrintInput();
#endif

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

	//runMopac(readInp_, allAtoms);

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
	h2o[0].atomlabel = "H";
	h2o[0].x = 0.757e0;
	h2o[0].y = 0.586e0;
	h2o[0].z = 0.0e0;
	h2o[0].atomlabel = "H";
	h2o[0].x = -0.757e0;
	h2o[0].y = 0.586e0;
	h2o[0].z = 0.0e0;
	string titleInfoh2o = "H2O monodentate 2";
	H2o_.setLigandCoordinates(h2o, titleInfoh2o);
}



