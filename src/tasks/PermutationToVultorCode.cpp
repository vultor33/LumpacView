#include "PermutationToVultorCode.h"

#include <iostream>

#include "CauchyIndex.h"

using namespace std;

PermutationToVultorCode::PermutationToVultorCode(){}
PermutationToVultorCode::~PermutationToVultorCode(){}

/*
INPUT INFORMATION NEEDED
* geo_code (consult Geometris.cpp)
* composition - e.x. Ma4(AB)
* permutation - e.x. 0 2 3 1
OUTPUT
* point group (string)

HOW TO RUN

1. Transform composition to accepted format:
   e.x. Maaabbc = Ma3b2c

2. Create a file in ./data with the name:
   "GEO-XX-" + your_composition + ".csv"
   e.x. GEO-10-Maaaaaaaaaa.csv

3. Put atom types in the first line. (start with 1)
   e.x. 2 3 3 1 1 

4. Put your permutation starting with 1,
   following the template:
   [DUM] [your_permutation_here] A {DUMMY}
   e.x. [DUM] [1 2 3 4 5 6 7 8 9 10] A {DUMMY}

5. Hard code this file for parameters

6. Run this task. All your permutations will be written in
   the .exe folder with the name input_name.
*/

void PermutationToVultorCode::run()
{
    int geo_code = 103; // PUT YOUR GEO CODE HERE
    string input_name = "GEO-10-Maaaaaaaaaa.csv"; // PUT YOUR FILE NAME HERE
    string input_path = "/home/fred/pesquisa/LumpacView/data/"; //PUT THE PATH OF YOUR DATA FOLDER HERE

	CauchyIndex ciSymmetry_(geo_code);
	ciSymmetry_.findAllSymmetryOperations(geo_code, input_name, input_path);
}
