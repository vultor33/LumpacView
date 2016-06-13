#ifndef ROOTMEANSQUAREDEVIATION_H
#define ROOTMEANSQUAREDEVIATION_H

#include <string>
#include <vector>

#include "Coordstructs.h"

class RootMeanSquareDeviation
{
public:
	RootMeanSquareDeviation();
	
	~RootMeanSquareDeviation();

	double rmsOverlay(std::string molName1, std::string molName2);

	double rmsOverlay(std::vector<CoordXYZ> & mol1, std::vector<CoordXYZ> & mol2);

	std::vector<CoordXYZ> readCoord(std::string fName);

private:
	std::vector<double> readPoint(std::string fName, int format = 1);
	std::vector<double> rotateToZ0(const std::vector<double> &point);
	std::vector<double> rotateToPlane(const std::vector<double> &point);
	std::vector<double> mirrorY(const std::vector<double> &point);
	void printXyz(std::string fName, const std::vector<double> &points);
	void printXyz(std::string fName, std::vector<CoordXYZ> &mol);
	void printXyzSuperpositions(std::string fName, std::vector<CoordXYZ> &mol1, std::vector<CoordXYZ> &mol2);
	double rms(std::vector<CoordXYZ> &mol1, std::vector<CoordXYZ> &mol2);
	void rotateMol(std::vector<CoordXYZ> &mol, double x, double y, double z, double angle);

	double oldRmsdToAddressPoints(std::string file1, std::string file2);
	double rmsOverlayThrouhgThreePoints(std::vector<CoordXYZ> & mol1, std::vector<CoordXYZ> & mol2);

};


#endif

/*
TESTE COM TRIPSINA
Molecula 1
24
t
N     2.311   -8.843   -9.584
H     3.263   -8.506   -9.584
C     1.628   -8.360   -8.401
H     2.197   -8.717   -7.543
C     0.201   -8.886   -8.338
O    -0.241   -9.367   -7.297
C     1.563   -6.836   -8.388
H     2.571   -6.429   -8.464
H     0.790   -6.429   -9.040
C     1.372   -5.358   -8.147
C     1.777   -4.353   -8.935
H     2.312   -4.446   -9.881
C     0.737   -4.746   -7.058
N     1.375   -3.136   -8.301
H     1.540   -2.206   -8.659
C     0.757   -3.422   -7.182
C     0.130   -5.329   -5.939
H     0.097   -6.410   -5.808
C     0.201   -2.534   -6.253
H     0.239   -1.454   -6.393
C    -0.436   -4.454   -4.990
H    -0.917   -4.873   -4.106
C    -0.401   -3.113   -5.142
H    -0.852   -2.474   -4.383


Molecula 2
24
t
N     7.387    1.335    0.000
H     8.339    1.672   -0.000
C     6.704    1.818    1.183
H     7.273    1.460    2.041
C     5.277    1.291    1.246
O     4.835    0.811    2.287
C     6.639    3.342    1.196
H     7.647    3.749    1.120
H     6.030    3.617    0.335
C     5.943    3.940    2.395
C     5.751    5.242    2.649
H     6.069    6.078    2.026
C     5.363    3.260    3.473
N     5.049    5.326    3.892
H     4.762    6.181    4.346
C     4.840    4.116    4.346
C     5.279    1.884    3.715
H     5.697    1.156    3.019
C     4.190    3.739    5.528
H     3.775    4.473    6.218
C     4.628    1.481    4.900
H     4.544    0.417    5.120
C     4.107    2.373    5.768
H     3.612    2.016    6.672

*/