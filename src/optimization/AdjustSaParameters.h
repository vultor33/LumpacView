#ifndef ADJUSTSAPARAMETERS_H
#define ADJUSTSAPARAMETERS_H


class AdjustSaParameters
{
public:
	AdjustSaParameters(
		double saTemperatureUpdate_in,
		double maxAlfaAngle_in,
		double maxBetaAngle_in,
		double saInitialTemperature_in,
		double saAcceptance_in);

	~AdjustSaParameters();

	bool takeParametersFromFile();

	double saTemperatureUpdate;
	double maxAlfaAngle;
	double maxBetaAngle;
	double saInitialTemperature;
	double saAcceptance;

private:

};

#endif
