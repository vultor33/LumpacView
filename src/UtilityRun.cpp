#include "UtilityRun.h"

#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "ChangeNames.h"
#include "CauchyIndex.h"
#include "IdentifyIsomers.h"
#include "Geometries.h"
#include "Coordstructs.h"
#include "IsomersToMol.h"
#include "ReadWriteFormats.h"


using namespace std;

UtilityRun::UtilityRun(){}

UtilityRun::~UtilityRun(){}

void UtilityRun::renameAtomTypes(string responseName)
{
	ifstream response_(responseName.c_str());
	string line;
	while (!response_.eof())
	{
		getline(response_, line);
		if (line == "")
			break;

		stringstream convert;
		convert << line;
		string comp;
		convert >> comp;
		string oldName = comp + "---atomTypes.txt";
		string newName = "final-" + comp + "-atomTypes";
		int result = rename(oldName.c_str(), newName.c_str());
		if (result == 0)
		{
			cout << comp << " ...  ...  ...  Ok" << endl;
		}
		else
		{
			cout << "Erro na hora de renomear" << endl;
			exit(1);
		}
	}
}

void UtilityRun::formatToSymmetryAndFiles(int geoCode)
{
	Geometries geo_;
	vector<CoordXYZ> mol0;
	double cutAngle;
	vector<int> reflec;
	geo_.selectGeometry(geoCode, mol0, cutAngle, reflec);
	string responseName = getResponseName(mol0.size());

	cout << "antes do response" << endl;
	
	findAllGroupPoint(geoCode);

	cout << "find group done" << endl;

	ChangeNames chnamessda;
	chnamessda.createNewCounting(geoCode, "", responseName);
	/*
    	IsomersToMol ismol_;
		ismol_.printAllMolFromSpecifiedGeometry(geoCode,
		"",
		responseName);
	*/

}

void UtilityRun::formatIsomersFiles()
{
	
	ChangeNames chNames_;
	string responseName = "response-combinations4.txt";
	/*
	chNames_.changeNameOfFiles(
		responseName,
		40,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#4\\T\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#4\\T\\");
	chNames_.changeNameOfFiles(
		responseName,
		41,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#4\\SP\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#4\\SP\\");

	responseName = "response-combinations5.txt";
	chNames_.changeNameOfFiles(
		responseName,
		50,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#5\\TBPY\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#5\\TBPY\\");
	chNames_.changeNameOfFiles(
		responseName,
		51,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#5\\SPY\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#5\\SPY\\");

	responseName = "response-combinations6.txt";
	chNames_.changeNameOfFiles(
		responseName,
		60,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#6\\OC\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#6\\OC\\");
	chNames_.changeNameOfFiles(
		responseName,
		61,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#6\\TPR\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#6\\TPR\\");

	responseName = "response-combinations7.txt";
	chNames_.changeNameOfFiles(
		responseName,
		70,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#7\\COC\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\COC\\");
	chNames_.changeNameOfFiles(
		responseName,
		71,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#7\\PBPY\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\PBPY\\");
	chNames_.changeNameOfFiles(
		responseName,
		72,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#7\\CTPR\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\CTPR\\");

	responseName = "response-combinations8.txt";
	chNames_.changeNameOfFiles(
		responseName,
		80,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#8\\SAPR\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\SAPR\\");
	*/

	responseName = "response-combinations8.txt";
	chNames_.changeNameOfFiles(
		responseName,
		81,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#8\\TDD\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\TDD\\");

	/*
	chNames_.changeNameOfFiles(
		responseName,
		82,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#8\\BTPR\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\BTPR\\");
	chNames_.changeNameOfFiles(
		responseName,
		83,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#8\\HBPY\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\HBPY\\");
	chNames_.changeNameOfFiles(
		responseName,
		84,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#8\\CU\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\CU\\");

	responseName = "response-combinations9.txt";
	chNames_.changeNameOfFiles(
		responseName,
		90,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#9\\TCTPR\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\TCTPR\\");
	chNames_.changeNameOfFiles(
		responseName,
		91,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#9\\CSAPR\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\CSAPR\\");

	responseName = "response-combinations9.txt";
	chNames_.changeNameOfFiles(
		responseName,
		92,
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\RAW\\#9\\MFF\\",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\MFF\\");
		*/

}

void UtilityRun::identifyOne()
{
	IdentifyIsomers ident_;
	Geometries geom_;
	vector<CoordXYZ> mol0;
	vector<int> reflection;
	double cutAngle;
	int size;

	size = 60;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma2b2(AA).csv",
		"",
		"ERUJIM-cores(2).xyz",
		getCountingPath(size));


}


void UtilityRun::identifyAll()
{
	IdentifyIsomers ident_;
	Geometries geom_;
	vector<CoordXYZ> mol0;
	vector<int> reflection;
	double cutAngle;
	int size;

	/*
	size = 40;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"T-4-Ma3b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\T-4\\Nd-QAYWOC\\",
		"QAYWOC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"T-4-Ma3b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\T-4\\Yb-HOJFAN\\",
		"HOJFAN.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"T-4-Ma3b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\T-4\\Yb-IBIFII\\",
		"IBIFII.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"T-4-Ma3b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\T-4\\Yb-IBIGAB\\",
		"IBIGAB.xyz",
		getCountingPath(size));


	size = 50;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"TBPY-5-Ma3b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TBPY-5\\La-LAPTEB\\",
		"LAPTEB10.xyz",
		getCountingPath(size));

	size = 51;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"SPY-5-Ma3b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SPY-5\\La-ZIDSOX\\",
		"ZIDSOX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SPY-5-Ma3b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SPY-5\\Lu-POGWEN\\",
		"POGWEN.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SPY-5-Ma3b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SPY-5\\Pr-POGWIR\\",
		"POGWIR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SPY-5-Ma3b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SPY-5\\Sm-NAFKIO\\",
		"NAFKIO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SPY-5-Ma3b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SPY-5\\Yb-JEMROI\\",
		"JEMROI.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SPY-5-Ma(AB)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SPY-5\\Yb-ZALFAT\\",
		"ZALFAT.xyz",
		getCountingPath(size));

	size = 60;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma3b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Dy-WAQZEU\\",
		"WAQZEU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Dy-ZAXSAS\\",
		"ZAXSAS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma3b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Er-ZADWUW\\",
		"ZADWUW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma4b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Eu-NURFIP\\",
		"NURFIP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Gd-WEWNOB\\",
		"WEWNOB.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma4b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Ho-NUYNOL\\",
		"NUYNOL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma4b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Nd-MINLIE\\",
		"MINLIE.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AA)2(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Pm-NUQYUT\\",
		"NUQYUT.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Pr-KUSXUR\\",
		"KUSXUR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Pr-ZAXRUL\\",
		"ZAXRUL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma2b2cd.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Sm-XUYPUC\\",
		"XUYPUC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma2b2(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Tb-BUJCAL\\",
		"BUJCAL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma5b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Tb-IXOBIG\\",
		"IXOBIG.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma6.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Tb-VAPTEL\\",
		"VAPTEL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-Ma6.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Tb-VAPTEL01\\",
		"VAPTEL01.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Yb-OFOPII\\",
		"OFOPII.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"OC-6-M(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\OC-6\\Yb-OFOPOO\\",
		"OFOPOO.xyz",
		getCountingPath(size));

	size = 61;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"TPR-6-M(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TPR-6\\Dy-USEPEO\\",
		"USEPEO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TPR-6-M(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TPR-6\\Dy-XAYRIZ\\",
		"XAYRIZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TPR-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TPR-6\\Er-TMHDER\\",
		"TMHDER.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TPR-6-M(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TPR-6\\Er-XOVHAS\\",
		"XOVHAS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TPR-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TPR-6\\Gd-YUWZOF\\",
		"YUWZOF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TPR-6-M(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TPR-6\\Yb-RENXIR\\",
		"RENXIR.xyz",
		getCountingPath(size));

	size = 70;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"COC-7-Ma(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\COC-7\\Er-XOYXIS\\",
		"XOYXIS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"COC-7-Ma(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\COC-7\\Eu-QHDOEU\\",
		"QHDOEU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"COC-7-Ma(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\COC-7\\Ho-PHPRHO10\\",
		"PHPRHO10.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"COC-7-Ma(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\COC-7\\Sm-XAXYAW\\",
		"XAXYAW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"COC-7-Ma(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\COC-7\\Tb-XAXXUP\\",
		"XAXXUP.xyz",
		getCountingPath(size));

	size = 71;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Dy-RABBEX\\",
		"RABBEX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Er-RIKTEK\\",
		"RIKTEK.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Er-RIKTEK01\\",
		"RIKTEK01.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma4b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Eu-VIGPAC\\",
		"VIGPAC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3b2(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Gd-GOZWEX\\",
		"GOZWEX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Ho-EWIPUV\\",
		"EWIPUV.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Nd-JALNAM\\",
		"JALNAM.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma4b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Nd-LUSNOD\\",
		"LUSNOD.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3b2c2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Nd-NATDET\\",
		"NATDET.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Nd-XIGKUT\\",
		"XIGKUT.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma4b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Pm-Pm_VIGPAC\\",
		"Pm_VIGPAC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Sm-AQEZIF\\",
		"AQEZIF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Sm-AQEZIF01\\",
		"AQEZIF01.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma4b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Sm-JALMIT\\",
		"JALMIT.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Sm-JALNEQ\\",
		"JALNEQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AB)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Sm-OFULAC\\",
		"OFULAC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Tb-EWIPOP\\",
		"EWIPOP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"PBPY-7-Ma(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\PBPY-7\\Yb-QAKXIJ\\",
		"QAKXIJ.xyz",
		getCountingPath(size));

	size = 72;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"CTPR-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CTPR-7\\Eu-JALNIU\\",
		"JALNIU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CTPR-7-Ma3(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CTPR-7\\Nd-SAXJIL01\\",
		"SAXJIL01.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CTPR-7-Ma(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CTPR-7\\Pr-ZAXSEW\\",
		"ZAXSEW.xyz",
		getCountingPath(size));

	size = 80;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ce-NOJTAH\\",
		"NOJTAH-PARTE-1-IDENTIFICACAO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ce-NOJTAH\\",
		"NOJTAH-PARTE-2-IDENTIFICACAO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ce-NOJTEL\\",
		"NOJTEL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ce-PEKWEH\\",
		"PEKWEH.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ce-PUTQAW\\",
		"PUTQAW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Dy-GAKYEW\\",
		"GAKYEW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Dy-PALBIN\\",
		"PALBIN.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Er-DIYNII\\",
		"DIYNII.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Er-GAKYOG\\",
		"GAKYOG.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Er-XEWWUR\\",
		"XEWWUR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2b2(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Er-ZUFSAU\\",
		"ZUFSAU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-ACPNEU\\",
		"ACPNEU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-JAXXOV\\",
		"JAXXOV.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-KELNOE\\",
		"KELNOE.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-LEJTAV\\",
		"LEJTAV.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-LELRUP\\",
		"LELRUP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-MIHNOG\\",
		"MIHNOG.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-PIEUAC01\\",
		"PIEUAC01.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-QALFOY\\",
		"QALFOY.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-TMHPEU10\\",
		"TMHPEU10.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-YEZFAK\\",
		"YEZFAK.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2b2(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-YODYIZ\\",
		"YODYIZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-YOJDIK\\",
		"YOJDIK.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2b2(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Eu-ZEXJUH\\",
		"ZEXJUH.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Gd-DIYNEE\\",
		"DIYNEE.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Gd-DUFBEL\\",
		"DUFBEL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Gd-GAKYAS\\",
		"GAKYAS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Gd-LEJVEB\\",
		"LEJVEB.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ho-FAGYOC\\",
		"FAGYOC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ho-GAKYIA\\",
		"GAKYIA.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Ho-SIFZIQ\\",
		"SIFZIQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma8.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Lu-FOPPOP\\",
		"FOPPOP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Lu-SUDDOK\\",
		"SUDDOK.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma3b(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Nd-DICNUZ\\",
		"DICNUZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Nd-LEJSUO\\",
		"LEJSUO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Mab(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Nd-XIPKIQ\\",
		"XIPKIQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Pr-CAZGUF\\",
		"CAZGUF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Pr-DIYMUT\\",
		"DIYMUT.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Pr-LEJSOI\\",
		"LEJSOI.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Pr-PELGOC\\",
		"PELGOC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Sm-CAZHAM\\",
		"CAZHAM.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Sm-FUJYEO\\",
		"FUJYEO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma3b(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Sm-WIVBUZ\\",
		"WIVBUZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma7b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Tb-FAGZAP\\",
		"FAGZAP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Tb-LEJTEZ\\",
		"LEJTEZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Tb-PEJZAF\\",
		"PEJZAF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Tb-QALFUE\\",
		"QALFUE.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma8.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Tb-SEGVEF\\",
		"SEGVEF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Tm-MIHPAU\\",
		"MIHPAU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Yb-DIYNOO\\",
		"DIYNOO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"SAPR-8-Ma6b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\SAPR-8\\Yb-GAKYUM\\",
		"GAKYUM.xyz",
		getCountingPath(size));

	size = 81;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma4b4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Ce-AFURUO\\",
		"AFURUO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Dy-XIVFUD\\",
		"XIVFUD.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma4(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Dy-YAVSOD\\",
		"YAVSOD.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Er-MAGDOP\\",
		"MAGDOP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma6(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Er-XEWVOK\\",
		"XEWVOK.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Er-YEGFEV\\",
		"YEGFEV.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Er-YEMSIT\\",
		"YEMSIT.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Eu-CIRKET\\",
		"CIRKET.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Eu-KAKPAN\\",
		"KAKPAN.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Eu-MASKAS\\",
		"MASKAS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Eu-QAKWUU\\",
		"QAKWUU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Eu-XIWTUS\\",
		"XIWTUS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Gd-FONMEA\\",
		"FONMEA.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Gd-NIVQEO\\",
		"NIVQEO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)3(AB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Gd-XIVFOX\\",
		"XIVFOX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Gd-ZIZNUR\\",
		"ZIZNUR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma6(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Ho-XEWVIE\\",
		"XEWVIE.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma4b2(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Ho-XORGEQ\\",
		"XORGEQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\La-QAKWEE\\",
		"QAKWEE.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Lu-XEPLUZ\\",
		"XEPLUZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma6(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Lu-XEWWAX\\",
		"XEWWAX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Nd-GOPRIM\\",
		"GOPRIM.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Nd-HEBCIA\\",
		"HEBCIA.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Nd-JOCDAH\\",
		"JOCDAH.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Pm-FUJYEO\\",
		"FUJYEO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma4(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Pm-GUPHUU\\",
		"GUPHUU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Sm-WOKHUZ\\",
		"WOKHUZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-M(AA)3(AB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Sm-XIVFIR\\",
		"XIVFIR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Tb-XARXET\\",
		"XARXET.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma4b2(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Tb-XORGAM\\",
		"XORGAM.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TDD-8-Ma6(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TDD-8\\Yb-XEWVUQ\\",
		"XEWVUQ.xyz",
		getCountingPath(size));

	size = 82;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma3b(AB)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Er-AKIYEY\\",
		"AKIYEY.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Er-DIJQAO\\",
		"DIJQAO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Er-TEJFEU\\",
		"TEJFEU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Er-ZOFJUZ\\",
		"ZOFJUZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Eu-NOJMOO01\\",
		"NOJMOO01.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma5b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Eu-PITCUQ\\",
		"PITCUQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Eu-TEJDES\\",
		"TEJDES.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma4(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Ho-BEYSAZ\\",
		"BEYSAZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Nd-GAGHED\\",
		"GAGHED.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Nd-TEJDIW\\",
		"TEJDIW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Pm-CAZHAM\\",
		"CAZHAM.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Sm-AXUROB\\",
		"AXUROB.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma4(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Sm-GUPHUU\\",
		"GUPHUU.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-M(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Sm-TEJDOC\\",
		"TEJDOC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma2(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Sm-TIMPUA\\",
		"TIMPUA.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"BTPR-8-Ma7b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\BTPR-8\\Tm-FAGYUI\\",
		"FAGYUI.xyz",
		getCountingPath(size));

	size = 83;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"HBPY-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\HBPY-8\\Er-DIJQIW\\",
		"DIJQIW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"HBPY-8-Ma2(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\HBPY-8\\Er-WEFVIM\\",
		"WEFVIM.xyz",
		getCountingPath(size));

	size = 84;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"CU-8-M(AA)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CU-8\\Sm-XEPLAF\\",
		"XEPLAF.xyz",
		getCountingPath(size));

	size = 90;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma9.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\a9\\Ce-CIBSAH\\",
		"CIBSAH.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Dy-CECLIF\\",
		"CECLIF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Dy-CECLIF10\\",
		"CECLIF10.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Er-KOZBUW\\",
		"KOZBUW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Er-VUSGUL\\",
		"VUSGUL.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma2b(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Eu-QIQHAZ\\",
		"QIQHAZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Eu-VUSGOF\\",
		"VUSGOF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\La-GOJQAX\\",
		"GOJQAX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma7b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\La-GULFOI\\",
		"GULFOI.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Nd-ANTNND10\\",
		"ANTNND10.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma7(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Nd-CAHJAX\\",
		"CAHJAX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Nd-TUPYOS\\",
		"TUPYOS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma8b.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Pr-FAGYIW\\",
		"FAGYIW.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma7(AA).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Pr-XOKYIF\\",
		"XOKYIF.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Sm-MOXJEO\\",
		"MOXJEO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Sm-YENHOO\\",
		"YENHOO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Tb-TOKVIY\\",
		"TOKVIY.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"TCTPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\TCTPR-9\\Tm-TUPYUY\\",
		"TUPYUY.xyz",
		getCountingPath(size));

	size = 91;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma(AB)4.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Eu-ZACXAC\\",
		"ZACXAC.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma7b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Gd-PEBDOP\\",
		"PEBDOP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma6b3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\La-XEMNUY\\",
		"XEMNUY.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma2b(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Lu-POHDIZ\\",
		"POHDIZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma3(AB)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Nd-ZANSIS\\",
		"ZANSIS.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Pm-QALFAK\\",
		"QALFAK.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma7b2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Pr-QOVXII\\",
		"QOVXII.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"CSAPR-9-Ma(AA)2(BB)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\CSAPR-9\\Yb-WATZEW\\",
		"WATZEW.xyz",
		getCountingPath(size));
	*/

	size = 92;
	geom_.selectGeometry(size, mol0, cutAngle, reflection);
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Dy-DIBTIR\\",
		"DIBTIR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Er-DIBTAJ\\",
		"DIBTAJ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma2b(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Er-DOGKEP\\",
		"DOGKEP.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Er-SEGVAB\\",
		"SEGVAB.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Er-TUMJEQ\\",
		"TUMJEQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma2b(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Er-VOSNOG\\",
		"VOSNOG.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma2b(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Eu-GEBYAN\\",
		"GEBYAN.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma2b(AA)2(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Eu-YUXREO\\",
		"YUXREO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma(AA)3(BB).csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Gd-JARBUZ\\",
		"JARBUZ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Gd-LANITB\\",
		"LANITB.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma2b(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Gd-NAVWIQ\\",
		"NAVWIQ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\La-SUXLIG\\",
		"SUXLIG.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma5(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\La-ZEJFOJ\\",
		"ZEJFOJ.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Lu-FEWKEX\\",
		"FEWKEX.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Lu-HELGUA\\",
		"HELGUA.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Nd-CUYLEO\\",
		"CUYLEO.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma4b(AA)2.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Nd-YENKOR\\",
		"YENKOR.xyz",
		getCountingPath(size));
	ident_.coordinatesToPermutation(
		mol0,
		"MFF-9-Ma3(AA)3.csv",
		"C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\IDENTIFICACAO-raio-x\\MFF-9\\Sm-QALFAK\\",
		"QALFAK.xyz",
		getCountingPath(size));


}

void UtilityRun::generateAllIsomersMol2Files()
{
	IsomersToMol ismol40_;
	string responseName = "response-combinations4.txt";
	ismol40_.printAllMolFromSpecifiedGeometry(40,
		getResultsPath(40),
		responseName);
	IsomersToMol ismol41_;
	responseName = "response-combinations4.txt";
	ismol41_.printAllMolFromSpecifiedGeometry(41,
		getResultsPath(41),
		responseName);

	IsomersToMol ismol50_;
	responseName = "response-combinations5.txt";
	ismol50_.printAllMolFromSpecifiedGeometry(50,
		getResultsPath(50),
		responseName);
	IsomersToMol ismol51_;
	responseName = "response-combinations5.txt";
	ismol51_.printAllMolFromSpecifiedGeometry(51,
		getResultsPath(51),
		responseName);

	IsomersToMol ismol60_;
	responseName = "response-combinations6.txt";
	ismol60_.printAllMolFromSpecifiedGeometry(60,
		getResultsPath(60),
		responseName);
	IsomersToMol ismol61_;
	responseName = "response-combinations6.txt";
	ismol61_.printAllMolFromSpecifiedGeometry(61,
		getResultsPath(61),
		responseName);

	IsomersToMol ismol70_;
	responseName = "response-combinations7.txt";
	ismol70_.printAllMolFromSpecifiedGeometry(70,
		getResultsPath(70),
		responseName);
	IsomersToMol ismol71_;
	responseName = "response-combinations7.txt";
	ismol71_.printAllMolFromSpecifiedGeometry(71,
		getResultsPath(71),
		responseName);
	IsomersToMol ismol72_;
	responseName = "response-combinations7.txt";
	ismol72_.printAllMolFromSpecifiedGeometry(72,
		getResultsPath(72),
		responseName);

	IsomersToMol ismol80_;
	responseName = "response-combinations8.txt";
	ismol80_.printAllMolFromSpecifiedGeometry(80,
		getResultsPath(80),
		responseName);
	IsomersToMol ismol81_;
	responseName = "response-combinations8.txt";
	ismol81_.printAllMolFromSpecifiedGeometry(81,
		getResultsPath(81),
		responseName);
	IsomersToMol ismol82_;
	responseName = "response-combinations8.txt";
	ismol82_.printAllMolFromSpecifiedGeometry(82,
		getResultsPath(82),
		responseName);
	IsomersToMol ismol83_;
	responseName = "response-combinations8.txt";
	ismol83_.printAllMolFromSpecifiedGeometry(83,
		getResultsPath(83),
		responseName);
	IsomersToMol ismol84_;
	responseName = "response-combinations8.txt";
	ismol84_.printAllMolFromSpecifiedGeometry(84,
		getResultsPath(84),
		responseName);

	IsomersToMol ismol90_;
	responseName = "response-combinations9.txt";
	ismol90_.printAllMolFromSpecifiedGeometry(90,
		getResultsPath(90),
		responseName);
	IsomersToMol ismol91_;
	responseName = "response-combinations9.txt";
	ismol91_.printAllMolFromSpecifiedGeometry(91,
		getResultsPath(91),
		responseName);
	IsomersToMol ismol92_;
	responseName = "response-combinations9.txt";
	ismol92_.printAllMolFromSpecifiedGeometry(92,
		getResultsPath(92),
		responseName);
}

void UtilityRun::findAllGroupPoint(int geoCode)
{
	cout << "ONLY ON WINDOWS  --- to linux, change: getResultsPath" << endl;

	string pathRead = getResultsPath(geoCode);
	Geometries geo_;
	vector<CoordXYZ> mol0;
	double cutAngle;
	vector<int> reflec;
	geo_.selectGeometry(geoCode, mol0, cutAngle, reflec);
	string responseName = getResponseName(mol0.size());

	ReadWriteFormats rwf_;
	string geomName = geo_.sizeToGeometryCode(geoCode);
	ifstream response_((pathRead + responseName).c_str());
	string line;

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

		cout << "loop find:  " << allIsomersCombinationFile << endl;

		CauchyIndex ciSymmetry_(geoCode);
		ciSymmetry_.findAllSymmetryOperations(
			geoCode,
			allIsomersCombinationFile,
			pathRead);
	}

}


string UtilityRun::getCountingPath(int geoCode)
{
	switch (geoCode)
	{
	case(40):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#4\\T\\T-4-counting.csv";
		break;

	case(41):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#4\\SP\\SP-4-counting.csv";
		break;

	case(50):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#5\\TBPY\\TBPY-5-counting.csv";
		break;

	case(51):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#5\\SPY\\SPY-5-counting.csv";
		break;

	case(60):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#6\\OC\\OC-6-counting.csv";
		break;

	case(61):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#6\\TPR\\TPR-6-counting.csv";
		break;

	case(70):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\COC\\COC-7-counting.csv";
		break;

	case(71):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\PBPY\\PBPY-7-counting.csv";
		break;

	case(72):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\CTPR\\CTPR-7-counting.csv";
		break;

	case(80):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\SAPR\\SAPR-8-counting.csv";
		break;

	case(81):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\TDD\\TDD-8-counting.csv";
		break;

	case(82):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\BTPR\\BTPR-8-counting.csv";
		break;

	case(83):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\HBPY\\HBPY-8-counting.csv";
		break;

	case(84):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\CU\\CU-8-counting.csv";
		break;

	case(90):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\TCTPR\\TCTPR-9-counting.csv";
		break;

	case(91):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\CSAPR\\CSAPR-9-counting.csv";
		break;

	case(92):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\MFF\\MFF-9-counting.csv";
		break;

	default:
		cout << "UtilityRun::getCountingPath geoCode not found" << endl;
		exit(1);
		break;
	}
}


std::string UtilityRun::getResultsPath(int geoCode)
{
	switch (geoCode)
	{
	case(40):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#4\\T\\";
		break;

	case(41):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#4\\SP\\";
		break;

	case(50):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#5\\TBPY\\";
		break;

	case(51):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#5\\SPY\\";
		break;

	case(60):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#6\\OC\\";
		break;

	case(61):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#6\\TPR\\";
		break;

	case(70):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\COC\\";
		break;

	case(71):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\PBPY\\";
		break;

	case(72):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#7\\CTPR\\";
		break;

	case(80):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\SAPR\\";
		break;

	case(81):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\TDD\\";
		break;

	case(82):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\BTPR\\";
		break;

	case(83):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\HBPY\\";
		break;

	case(84):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#8\\CU\\";
		break;

	case(90):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\TCTPR\\";
		break;

	case(91):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\CSAPR\\";
		break;

	case(92):
		return "C:\\Users\\basta\\Documents\\DOUTORADO\\!Trabalhos-Paralelo\\!QUALIFICACAO\\lumpac-view\\!!RESULTADOS\\FORMATADO\\#9\\MFF\\";
		break;

	default:
		cout << "UtilityRun::getCountingPath geoCode not found" << endl;
		exit(1);
		break;
	}
}


std::string UtilityRun::getResultsPathLinux(int geoCode)
{
	switch (geoCode)
	{
	case(40):
		return "//home//vultor//lumpacview//FORMATADO//#4//T//";
		break;

	case(41):
		return "//home//vultor//lumpacview//FORMATADO//#4//SP//";
		break;

	case(50):
		return "//home//vultor//lumpacview//FORMATADO//#5//TBPY//";
		break;

	case(51):
		return "//home//vultor//lumpacview//FORMATADO//#5//SPY//";
		break;

	case(60):
		return "//home//vultor//lumpacview//FORMATADO//#6//OC//";
		break;

	case(61):
		return "//home//vultor//lumpacview//FORMATADO//#6//TPR//";
		break;

	case(70):
		return "//home//vultor//lumpacview//FORMATADO//#7//COC//";
		break;

	case(71):
		return "//home//vultor//lumpacview//FORMATADO//#7//PBPY//";
		break;

	case(72):
		return "//home//vultor//lumpacview//FORMATADO//#7//CTPR//";
		break;

	case(80):
		return "//home//vultor//lumpacview//FORMATADO//#8//SAPR//";
		break;

	case(81):
		return "//home//vultor//lumpacview//FORMATADO//#8//TDD//";
		break;

	case(82):
		return "//home//vultor//lumpacview//FORMATADO//#8//BTPR//";
		break;

	case(83):
		return "//home//vultor//lumpacview//FORMATADO//#8//HBPY//";
		break;

	case(84):
		return "//home//vultor//lumpacview//FORMATADO//#8//CU//";
		break;

	case(90):
		return "//home//vultor//lumpacview//FORMATADO//#9//TCTPR//";
		break;

	case(91):
		return "//home//vultor//lumpacview//FORMATADO//#9//CSAPR//";
		break;

	case(92):
		return "//home//vultor//lumpacview//FORMATADO//#9//MFF//";
		break;

	default:
		cout << "UtilityRun::getCountingPath geoCode not found" << endl;
		exit(1);
		break;
	}
}


string UtilityRun::getResponseName(int size)
{
	switch (size)
	{
	case 4:
		return "response-combinations4.txt";
		break;

	case 5:
		return "response-combinations5.txt";
		break;

	case 6:
		return "response-combinations6.txt";
		break;

	case 7:
		return "response-combinations7.txt";
		break;

	case 8:
		return "response-combinations8.txt";
		break;

	case 9:
		return "response-combinations9.txt";
		break;

	case 10:
		return "response-combinations10.txt";
		break;

	case 11:
		return "response-combinations11.txt";
		break;

	case 12:
		return "response-combinations12.txt";
		break;

	default:
		cout << "ERROR ON UtilityRun::responseName - size not found" << endl;
		exit(1);
	}
}
