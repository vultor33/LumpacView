#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "MyExceptions.h"
#include "Coordstructs.h"
#include "ReadInput.h"
#include "ComplexCreator.h"
#include "ControlMopac.h"
#include "AuxMath.h"
#include "PointAnalysis.h"
#include "AdjustSaParameters.h"

using namespace std;

int main()
{
	// tinicial, saUpdate, maxAlfa, maxBeta,
	AdjustSaParameters saParameters_(
		0.13113111182697784,
		1.8245257230764378,
		2.0524859126887280,
		503.81980808157266,
		0.47976429810660104);

#ifdef _FITSA
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
#endif

#ifdef _DEBUG
		PointAnalysis p;
#endif

		ReadInput readInp_;
		try {
			readInp_.readLumpacViewInput();
		}
		catch (MyExceptions& e) {
			cout << e.what() << endl;
			return 1;
		}
		catch (...) {
			cout << "unknown error on input - check it or contact developers" << endl;
			return 1;
		}

#ifdef _DEBUG
		readInp_.rePrintInput();
#endif

		bool sucess;
		bool terminateExecution = false;
		for (size_t i = 0; i < readInp_.allLigands.size(); i++) {
			sucess = readInp_.allLigands[i].initializeLigand();
			if (!sucess) {
				cout << "Can't build ligand:  " << i << " check input" << endl;
				terminateExecution = true;
			}
		}
		if (terminateExecution) return 1;

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

		sucess = cpCreator.start();
		if (sucess) cout << "ligantes iniciados com sucesso" << endl;
		else {
			cout << "um problema aconteceu nos ligantes" << endl;
			return 1;
		}

		vector<CoordXYZ> allAtoms = cpCreator.simulatedAnnealing();

#ifdef _FITSA
		remove("fitness.txt");
		ofstream fit_("fitness.txt");
		fit_ << cpCreator.finalI << endl;
		fit_.close();
	}
	return 0;
#else
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
	sucess = controlMop.optimize(allAtoms);
	if (sucess) cout << "a estrutura otimizou com sucesso - tenha um bom dia" << endl;
	else {
		cout << "ocorreu algum problema na otimizacao do mopac, tente outra vez" << endl;
	}

	return 0;
#endif
}



/*
ERROR HANDLING
- Quando o destructor nao e chamado:
- Destructor of the class is not called if exception is thrown in its constructor.
- Exception is automatically re-thrown if caught in construction initialization list catch block.
*/


/*
to string
std::string to_string(const int& n);
string ReadInput::to_string( const int& n )
{
std::ostringstream stm ;
stm << n ;
return stm.str() ;
}





for (int i = 0; i < (int)ri_.numberLigands.size(); i++)
{
for (int j = 0; j < ri_.numberLigands[i]; j++)
{
Ligand auxOptLigand; //com 8 coordXYZ
for (int l = 0; l < (int)initialConfiguration[m].coord.size(); l++) //ligand i -> 8 atoms
{
m++;
CoordXYZ auxCoord;
auxCoord = optimizedMol_.coord[k];
k++;
auxCoord.x -= centerX;
auxCoord.y -= centerY;
auxCoord.z -= centerZ;
auxOptLigand.coord.push_back(auxCoord);
}
allOptimized.push_back(auxOptLigand); // add ligand i
}
}











for (int i = 0; i < 50; i++)
{

	Ligand h2o = ri_.readConfigurations("h2o-ligante");
	Ligand h2o2 = h2o;
	Ligand h2o3 = h2o;
	Ligand h2o4 = h2o;
	Ligand h2o5 = h2o;
	Ligand h2o6 = h2o;
	Ligand h2o7 = h2o;
	Ligand h2o8 = h2o;

	PositionChange pos_;
	pos_.initialPosition(h2o);
	pos_.initialPosition(h2o2);
	pos_.initialPosition(h2o3);
	pos_.initialPosition(h2o4);
	pos_.initialPosition(h2o5);
	pos_.initialPosition(h2o6);
	pos_.initialPosition(h2o7);
	pos_.initialPosition(h2o8);

	vector<Ligand> all(8);
	all[0] = h2o;
	all[1] = h2o2;
	all[2] = h2o3;
	all[3] = h2o4;
	all[4] = h2o5;
	all[5] = h2o6;
	all[6] = h2o7;
	all[7] = h2o8;

	pr_.printAllMol(all);
}
*/



