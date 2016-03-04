#include "AdjustSaParameters.h"

#include <fstream>

#include "AuxMath.h"

using namespace std;

AdjustSaParameters::AdjustSaParameters(
	double saTemperatureUpdate_in,
	double maxAlfaAngle_in,
	double maxBetaAngle_in,
	double saInitialTemperature_in,
	double saAcceptance_in)
{
	AuxMath auxMath_;
	saTemperatureUpdate = saTemperatureUpdate_in;
	maxAlfaAngle = maxAlfaAngle_in * (auxMath_._pi / 180.0e0);
	maxBetaAngle = maxBetaAngle_in * (auxMath_._pi / 180.0e0);
	saInitialTemperature = saInitialTemperature_in;
	saAcceptance = saAcceptance_in;
}

AdjustSaParameters::~AdjustSaParameters(){}

bool AdjustSaParameters::takeParametersFromFile()
{
	ifstream point_("point.txt");
	AuxMath auxMath_;
	int n;
	point_ >> n;
	point_ >> saTemperatureUpdate;
	if ((saTemperatureUpdate < 0) || (saTemperatureUpdate > 1))
		return false;

	point_ >> maxAlfaAngle;
	if ((maxAlfaAngle < 0) || (maxAlfaAngle > 90))
		return false;
	maxAlfaAngle *= (auxMath_._pi / 180.0e0);

	point_ >> maxBetaAngle;
	if ((maxBetaAngle < 0) || (maxBetaAngle > 90))
		return false;
	maxBetaAngle *= (auxMath_._pi / 180.0e0);

	point_ >> saInitialTemperature;
	if ((saInitialTemperature < 0) || (saInitialTemperature > 1000000))
		return false;

	point_ >> saAcceptance;
	if ((saAcceptance < 0) || (saAcceptance > 0.5e0))
		return false;

	point_.close();
	return true;

}