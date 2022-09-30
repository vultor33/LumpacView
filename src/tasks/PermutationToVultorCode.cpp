#include "PermutationToVultorCode.h"

#include <iostream>

#include "CauchyIndex.h"

using namespace std;

PermutationToVultorCode::PermutationToVultorCode(){}
PermutationToVultorCode::~PermutationToVultorCode(){}

void PermutationToVultorCode::run()
{
    std::cout << "Dentro do task" << std::endl;
    int geo_code = 101;
    string input_name = "GEO-10-Maaaaaaaaaa.csv";
    string input_path = "/home/fred/pesquisa/LumpacView/data/";

	CauchyIndex ciSymmetry_(geo_code);
	ciSymmetry_.findAllSymmetryOperations(geo_code, input_name, input_path);
}
