#ifndef AUXMATHGA_H
#define AUXMATHGA_H

#include <vector>
#include <fstream>

namespace zerg{
class AuxMathGa
{
public:
	//random numbers
	static void set_seed(int seed);
	static double randomNumber(double fMin, double fMax);
	static int randomNumber(int fMin, int fMax);
	static std::vector<double> unitarySphericalVector();

	//ordering
	static std::vector<int> vector_ordering(std::vector<double> &vetor_entrada); //organiza a entrada e solta um vetor de organizacao
	static void vector_ordering_with_instructions(std::vector<std::vector<double> > &vetor_entrada, const std::vector<int> &vetor_organiza); //organiza uma matriz com a instrucao obtida no anterior.
	static void vector_ordering_with_instructions(std::vector<int> &vetor_entrada, const std::vector<int> &vetor_organiza); //organiza uma matriz com a instrucao obtida no anterior.

	//linear algebra

};
}

#endif
