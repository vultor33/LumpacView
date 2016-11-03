#define useDlib

#include "DiagonalizeDlib.h"

#include <vector>



#ifdef useDlib
#include <dlib/matrix.h>
#include <dlib/matrix/matrix_eigenvalue.h>
using namespace dlib;
#endif

using namespace std;

DiagonalizeDlib::~DiagonalizeDlib(){}

DiagonalizeDlib::DiagonalizeDlib(std::vector< std::vector<double> > &entryMatrix)
{
#ifdef useDlib
	int size = entryMatrix.size();
	matrix<double> m;
	matrix<double> val;
	matrix<double> vec;
	m.set_size(size,size);
	vec.set_size(size,size);
	val.set_size(size,1);

	eigenvaluesDlib.resize(size);
	eigenvectorsDlib.resize(size);
	for(int i=0;i<size;i++)
	{
		eigenvectorsDlib[i].resize(size);
		for(int j=i; j<size; j++)
		{
			m(i,j) = entryMatrix[i][j];
			if(i!=j)
				m(j,i) = m(i,j);
		}
	}

	eigenvalue_decomposition< matrix<double> > eig_(m);
	val = eig_.get_real_eigenvalues();
	vec = eig_.get_pseudo_v();

	for(int i=0;i<size;i++)
	{
		eigenvaluesDlib[i] = val(i,0);		
		for(int j=0; j<size; j++)
		{
			eigenvectorsDlib[i][j] = vec(i,j);
		}
	}
#endif

}
