#include "AuxMath.h"

#include <cmath>
#include <stdlib.h>
#include <vector>

/* fredoldnormal - normalVectorFrom3Points - v1x = x1 - x2;	v1y = y1 - y2;	v1z = z1 - z2;	v2x = x1 - x3;	v2y = y1 - y3;	v2z = z1 - z3; e no rotateOnX1 o angle era + */

using namespace std;

AuxMath::~AuxMath(){}

AuxMath::AuxMath()
{
	_pi = 3.14159265358979323846;
}

vector<double> AuxMath::matrixXVector(vector< vector<double> > &M, double x, double y, double z)
{
	int size = M.size();
	vector<double> Mv(size);
	for (int i = 0; i < size; i++)
	{
		Mv[i] = M[i][0] * x + M[i][1] * y + M[i][2] * z;
	}
	return Mv;
}

double AuxMath::fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

int AuxMath::randomNumber(int fMin, int fMax)
{
	return fMin + (rand() % (int)(fMax - fMin + 1));
}

double AuxMath::norm(const double x, const double y, const double z)
{
	return sqrt(x * x + y * y + z * z);
}

void AuxMath::normalize(std::vector<double> &vector)
{
	double length = norm(vector[0], vector[1], vector[2]);
	vector[0] /= length;
	vector[1] /= length;
	vector[2] /= length;
}

vector<double> AuxMath::vectorProduct(double x1, double y1, double z1, double x2, double y2, double z2)
{
	vector<double> prodvec(3);
	prodvec[0] = y1*z2 - y2*z1;
	prodvec[1] = z1*x2 - z2*x1;
	prodvec[2] = x1*y2 - x2*y1;

	return prodvec;
}

vector< vector<double> > AuxMath::rotationMatrix(double x, double y, double z, double teta)
{
	vector< vector<double> > rot(3);
	rot[0].resize(3);
	rot[1].resize(3);
	rot[2].resize(3);
	double ct = cos(teta);
	double st = sin(teta);
	rot[0][0] = x*x*(1.0e0 - ct) + ct;
	rot[0][1] = x*y*(1.0e0 - ct) - z*st;
	rot[0][2] = x*z*(1.0e0 - ct) + y*st;
	rot[1][0] = x*y*(1.0e0 - ct) + z*st;
	rot[1][1] = y*y*(1.0e0 - ct) + ct;
	rot[1][2] = y*z*(1.0e0 - ct) - x*st;
	rot[2][0] = x*z*(1.0e0 - ct) - y*st;
	rot[2][1] = y*z*(1.0e0 - ct) + x*st;
	rot[2][2] = z*z*(1.0e0 - ct) + ct;
	return rot;
}

std::vector<double> AuxMath::normalVectorFrom3Points(
	double x1, double y1, double z1,
	double x2, double y2, double z2,
	double x3, double y3, double z3)
{
	double v1x, v1y, v1z, v2x, v2y, v2z;
	v1x = -x2 + x1;
	v1y = -y2 + y1;
	v1z = -z2 + z1;
	v2x = -x2 + x3;
	v2y = -y2 + y3;
	v2z = -z2 + z3;

	vector<double> n = vectorProduct(v1x, v1y, v1z, v2x, v2y, v2z);
	if ((abs(n[0]) < 1.0e-12) && (abs(n[1]) < 1.0e-12) && (abs(n[2]) < 1.0e-12))
		return n;

	normalize(n);
	return n;
}

std::vector<double> AuxMath::triangleCentroid(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
	vector<double> centroid(3);
	double terc = 1.0e0 / 3.0e0;
	centroid[0] = terc * (x1 + x2 + x3);
	centroid[1] = terc * (y1 + y2 + y3);
	centroid[2] = terc * (z1 + z2 + z3);

	return centroid;
}

double AuxMath::angleFrom3Points(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
{
	vector<double> v(3);
	vector<double> u(3);

	v[0] = -x2 + x1;
	v[1] = -y2 + y1;
	v[2] = -z2 + z1;
	u[0] = -x2 + x3;
	u[1] = -y2 + y3;
	u[2] = -z2 + z3;

	normalize(v);
	normalize(u);
	double prodInt = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
	if (prodInt < -1.0e0)
		prodInt = -1.0e0;
	else if (prodInt > 1.0e0)
		prodInt = 1.0e0;

	return acos(prodInt);
}

double AuxMath::escalarProduct(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return x1*x2 + y1*y2 + z1*z2;
}


vector<int> AuxMath::vector_ordering(vector<double> &vetor_entrada)
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
			if (aux_menor < vetor_entrada[j])
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

// USING
//vector<int> instructions = auxMath_.vector_ordering(vectorExemple);
//auxMath_.vector_ordering_with_instructions(otherVector, instructions);
vector<int> AuxMath::vector_ordering(vector<int> &vetor_entrada)
{
	int tamanho_vetor = vetor_entrada.size();
	vector<int> operacoes_organizam;

	int aux_menor;
	int aux_pos_menor;
	for (int i = 0; i < (tamanho_vetor - 1); i++)
	{
		aux_menor = vetor_entrada[i];
		aux_pos_menor = i;
		for (int j = i + 1; j < tamanho_vetor; j++)
		{
			if (aux_menor < vetor_entrada[j])
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
void AuxMath::vector_ordering_with_instructions(vector< vector<double> > &vetor_entrada, const vector<int> &vetor_organiza)
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

void AuxMath::vector_ordering_with_instructions(vector<int> &vetor_entrada, const vector<int> &vetor_organiza)
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

void AuxMath::vector_ordering_with_instructions(vector<bool> &vetor_entrada, const vector<int> &vetor_organiza)
{
	int tamanho_vetor_organiza = vetor_organiza.size();
	int tamanho_vetor_entrada = vetor_entrada.size();
	int aux_pos1, aux_pos2;
	bool aux_primeiro;

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


