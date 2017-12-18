#include "AuxMathGa.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace zerg;

namespace zerg{
void AuxMathGa::set_seed(int seed)
{
	srand(seed);
}

double AuxMathGa::randomNumber(double fMin, double fMax)
{
	double f = ((double)rand() / (double)(RAND_MAX));
  	return fMin + f * (fMax - fMin);
}

int AuxMathGa::randomNumber(int fMin, int fMax)
{
	return fMin + (rand() % (int)(fMax - fMin + 1));
}

vector<double> AuxMathGa::unitarySphericalVector()
{
	vector<double> unit(3);
	double pi = 3.1415926535897932384626433832795e0;
	double fi, teta;
	fi = 2.0e0 * pi * AuxMathGa::randomNumber(0.0e0, 1.0e0);
	teta = acos(2.0e0 * AuxMathGa::randomNumber(0.0e0, 1.0e0) - 1.0e0);
	unit[0] = sin(teta) * cos(fi);
	unit[1] = sin(teta) * sin(fi);
	unit[2] = cos(teta);
	return unit;
}

vector<int> AuxMathGa::vector_ordering(vector<double> &vetor_entrada)
{
	int tamanho_vetor = vetor_entrada.size();
	vector<int> operacoes_organizam;

	double aux_menor;
	int aux_pos_menor;
	for (int i = 0; i < (tamanho_vetor - 1); i++)
	{

		aux_menor = vetor_entrada[i];
		aux_pos_menor = i;
		for (int j = i + 1; j < tamanho_vetor; j++)
		{
			if (aux_menor > vetor_entrada[j])
			{
				aux_menor = vetor_entrada[j];
				aux_pos_menor = j;
			}
		}

		// Se nao, nem precisa anotar a troca.
		if (aux_pos_menor != i)
		{
			vetor_entrada[aux_pos_menor] = vetor_entrada[i];
			vetor_entrada[i] = aux_menor;
			operacoes_organizam.push_back(i);
			operacoes_organizam.push_back(aux_pos_menor);
		}

	}

	// Nessa saida esta gravado todas as informacoes que foram feitas para organizar
	// o vetor;
	return operacoes_organizam;
}

// Rotina que organiza os autovetores em ordem crescente de autovalores.
void AuxMathGa::vector_ordering_with_instructions(vector<vector<double> > &vetor_entrada, const vector<int> &vetor_organiza)
{
	int tamanho_vetor_organiza = vetor_organiza.size();
	int tamanho_vetor_entrada = vetor_entrada.size();
	int aux_pos1, aux_pos2;
	double aux_primeiro;

	for (int i = 0; i<tamanho_vetor_organiza; i += 2)
	{
		// o primeiro e igual ao segundo e o segundo e igual ao primeiro.
		aux_pos1 = vetor_organiza[i];
		aux_pos2 = vetor_organiza[i + 1];
		for (int j = 0; j < tamanho_vetor_entrada; j++)
		{
			aux_primeiro = vetor_entrada[j][aux_pos1];
			vetor_entrada[j][aux_pos1] = vetor_entrada[j][aux_pos2];
			vetor_entrada[j][aux_pos2] = aux_primeiro;
		}
	}

}

void AuxMathGa::vector_ordering_with_instructions(vector<int> &vetor_entrada, const vector<int> &vetor_organiza)
{
	int tamanho_vetor_organiza = vetor_organiza.size();
	int tamanho_vetor_entrada = vetor_entrada.size();
	int aux_pos1, aux_pos2;
	int aux_primeiro;

	for (int i = 0; i<tamanho_vetor_organiza; i += 2)
	{
		// o primeiro e igual ao segundo e o segundo e igual ao primeiro.
		aux_pos1 = vetor_organiza[i];
		aux_pos2 = vetor_organiza[i + 1];
		aux_primeiro = vetor_entrada[aux_pos1];
		vetor_entrada[aux_pos1] = vetor_entrada[aux_pos2];
		vetor_entrada[aux_pos2] = aux_primeiro;
	}
}

}
