/*
 *  *******************************************************************
 *
 *  rmsd.h 
 *  (c) 2005 Bosco K Ho
 * 
 *  Implementation of the Kabsch algorithm to find the RMSD, and 
 *  the least-squares rotation matrix for a superposition between 
 *  two sets of vectors.
 *
 *  This implementation is completely self-contained. No other dependencies.
 *
 *  **************************************************************************
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation; either version 2.1 of the License, or (at
 *  your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,  but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details. 
 *  
 *  You should have received a copy of the GNU Lesser General Public License 
 *  along with this program; if not, write to the Free Software Foundation, 
 *  Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 *  **************************************************************************
 * 
 */

#ifndef KABSCHRMSD_H
#define KABSCHRMSD_H

#include <vector>
#include <string>

#include "Coordstructs.h"

class KabschRmsd
{
public:
	KabschRmsd();
	
	~KabschRmsd();

	double rmsOverlay(std::string molName1, std::string molName2);

	double rmsOverlay(std::vector<CoordXYZ> & mol1, std::vector<CoordXYZ> & mol2);

	std::vector<CoordXYZ> readCoord(std::string fName);


private:
	double PI;

	/*
	* fast_rmsd()
	*
	* Fast calculation of rmsd w/o calculating a rotation matrix,
	* adapted from the BTK by Chris Saunders 11/2002.
	*/
	void fast_rmsd(double ref_xlist[][3],
		double mov_xlist[][3],
		int n_list,
		double* rmsd);


	/*
	 * calculate_rotation_rmsd()
	 *
	 *   given two lists of x,y,z coordinates, constructs
	 *    - mov_com: the centre of mass of the mov list
	 *    - mov_to_ref: vector between the com of mov and ref
	 *    - U: the rotation matrix for least-squares, usage of
	 *         of the matrix U[3][3] is
	 *           for (i=0; i<3; i++)
	 *           {
	 *             rotated_v[i] = 0.0;
	 *             for (j=0; j<3; j++)
	 *               rotated_v[i] += U[i][j] * v[j];
	 *           }
	 *    - rmsd: measures similarity between the vectors
	 */
	void calculate_rotation_rmsd(double ref_xlist[][3],
		double mov_xlist[][3],
		int n_list,
		double mov_com[3],
		double mov_to_ref[3],
		double U[3][3],
		double* rmsd);

	int calculate_rotation_matrix(double R[3][3],
		double U[3][3],
		double E0,
		double* residual);

	void normalize(double a[3]);
	
	double dot(double a[3], double b[3]);

	void cross(double a[3], double b[3], double c[3]);

	/*
	* setup_rotation()
	*
	*      given two lists of x,y,z coordinates, constructs
	* the correlation R matrix and the E value needed to calculate the
	* least-squares rotation matrix.
	*/
	void setup_rotation(double ref_xlist[][3],
		double mov_xlist[][3],
		int n_list,
		double mov_com[3],
		double mov_to_ref[3],
		double R[3][3],
		double* E0);

	int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot);

	int diagonalize_symmetric(double matrix[3][3],
		double eigen_vec[3][3],
		double eigenval[3]);


};

#endif


/*
EXEMPLO
vector<CoordXYZ> mol1(3);
vector<CoordXYZ> mol2(3);
mol1[0].x = 0;
mol1[0].y = 0;
mol1[0].z = 0;
mol1[1].x = 1;
mol1[1].y = 0;
mol1[1].z = 0;
mol1[2].x = 0;
mol1[2].y = 1;
mol1[2].z = 0;
mol2[0].x = 0;
mol2[0].y = 0;
mol2[0].z = 0;
mol2[1].x = 1;
mol2[1].y = 0;
mol2[1].z = 0;
mol2[2].x = 0;
mol2[2].y = 1;
mol2[2].z = 0;

KabschRmsd krmsd_;
cout << "rmsd:  " << krmsd_.rmsOverlay(mol1, mol2) << endl;


*/
