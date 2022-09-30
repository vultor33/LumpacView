#include "GroupPointIdentify.h"

#include <string>
#include <sstream>
#include <iostream>

#include "CauchyIndex.h"
#include "Geometries.h"
#include "ReadWriteFormats.h"

using namespace std;

GroupPointIdentify::GroupPointIdentify(){}

GroupPointIdentify::~GroupPointIdentify(){}

string GroupPointIdentify::findGroupPoint(vector<string> & allSymmetries)
{
	stringstream group;
	int princAxis = findPrincipalAxis(allSymmetries);
	if (princAxis == -1)
	{
		if (findPlane(allSymmetries))
			group << "Cs";
		else if (findInversion(allSymmetries))
			group << "Ci";
		else
			group << "C1";

		return group.str();
	}
	else if (rotationNumber(allSymmetries[princAxis]) > 2)
	{
		int otherCn = findOtherCnAxis(allSymmetries, princAxis);
		if (otherCn != -1)
		{
			if (rotationNumber(allSymmetries[otherCn]) > 2)
			{
				if (!findInversion(allSymmetries))
					group << "Td";
				else if (rotationNumber(allSymmetries[otherCn]) == 5)
					group << "Ih";
				else
					group << "Oh";

				return group.str();
			}
		}
	}

	if (findnC2Perpendicular(allSymmetries, princAxis))
	{
		if (findPlanePerpendicular(allSymmetries, princAxis))
			group << "D" << rotationNumber(allSymmetries[princAxis]) << "h";
		else if (findnPlanes(allSymmetries, princAxis))
			group << "D" << rotationNumber(allSymmetries[princAxis]) << "d";
		else
			group << "D" << rotationNumber(allSymmetries[princAxis]);

		return group.str();
	}
	else
	{
		if (findPlanePerpendicular(allSymmetries, princAxis))
			group << "C" << rotationNumber(allSymmetries[princAxis]) << "h";
		else if (findnPlanes(allSymmetries, princAxis))
			group << "C" << rotationNumber(allSymmetries[princAxis]) << "v";
		else if (find2nSimproper(allSymmetries, princAxis))
			group << "S" << 2 * rotationNumber(allSymmetries[princAxis]);
		else
			group << "C" << rotationNumber(allSymmetries[princAxis]);

		return group.str();
	}
}


int GroupPointIdentify::findPrincipalAxis(vector<string> & allSymmetries)
{
	int maxRot = 0;
	int maxPos = -1;
	for (int i = 0; i < (int)allSymmetries.size(); i++)
	{
		int rotI = rotationNumber(allSymmetries[i]);
		if (rotI > maxRot)
		{
			maxRot = rotI;
			maxPos = i;
		}
	}
	return maxPos;
}

bool GroupPointIdentify::findPlane(vector<string> & allSymmetries)
{
	for (size_t i = 0; i < allSymmetries.size(); i++)
	{
		if (allSymmetries[i][0] == 'P')
			return true;
	}
	return false;
}

bool GroupPointIdentify::findInversion(vector<string> & allSymmetries)
{
	for (size_t i = 0; i < allSymmetries.size(); i++)
	{
		if (allSymmetries[i][0] == 'I')
			return true;
	}
	return false;
}

int GroupPointIdentify::findOtherCnAxis(
	vector<string> & allSymmetries,
	int princAxis)
{
	string princAxisName = rotationCode(allSymmetries[princAxis]);
	int maxRot = 0;
	int maxPos = -1;
	for (int i = 0; i < (int)allSymmetries.size(); i++)
	{
		int rotI = rotationNumber(allSymmetries[i]);
		string axisI = rotationCode(allSymmetries[i]);
		if (axisI != princAxisName)
		{
			if (axisI[0] == 'C')
			{
				if ((axisI[3] != princAxisName[3]) ||
					(axisI[4] != princAxisName[4]))
				{
					if (rotI > maxRot)
					{
						maxRot = rotI;
						maxPos = i;
					}
				}
			}
		}
	}
	return maxPos;
}


bool GroupPointIdentify::findnC2Perpendicular(
	std::vector<std::string> & allSymmetries,
	int princAxis)
{
	int count = 0;
	for (size_t i = 0; i < allSymmetries.size(); i++)
	{
		if (i == princAxis)
			continue;
		if (rotationNumber(allSymmetries[i]) == 2)
		{
			if (findCodeOnParentheses(allSymmetries[i], allSymmetries[princAxis]))
				count++;
		}
	}

	return count >= rotationNumber(allSymmetries[princAxis]);
}

bool GroupPointIdentify::find2nSimproper(
	std::vector<std::string> & allSymmetries,
	int princAxis)
{
	int n = rotationNumber(allSymmetries[princAxis]);
	for (size_t i = 0; i < allSymmetries.size(); i++)
	{
		if (sImproperRotationNumber(allSymmetries[i]) == 2 * n)
			return true;
	}
	return false;
}

bool GroupPointIdentify::findPlanePerpendicular(
	std::vector<std::string> & allSymmetries,
	int princAxis)
{
	for (size_t i = 0; i < allSymmetries.size(); i++)
	{
		if (i == princAxis)
			continue;
		if (allSymmetries[i][0] == 'P')
		{
			if (findCodeOnParentheses(allSymmetries[i], allSymmetries[princAxis]))
				return true;
		}
	}
	return false;
}

bool GroupPointIdentify::findnPlanes(
	std::vector<std::string> & allSymmetries,
	int princAxis)
{
	int n = rotationNumber(allSymmetries[princAxis]);
	int count = 0;
	for (size_t i = 0; i < allSymmetries.size(); i++)
	{
		if (i == princAxis)
			continue;
		if (allSymmetries[i][0] == 'P')
		{
			count++;
		}
	}
	return count >= n;
}


string GroupPointIdentify::codeOnly(string wholeCode)
{
	stringstream convert;
	convert << wholeCode;
	string first;
	convert >> first;
	return first;
}


int GroupPointIdentify::rotationNumber(string code)
{
	if (code[0] == 'C')
	{
		stringstream convert;
		convert << code[1];
		int rotNumber;
		convert >> rotNumber;
		return rotNumber;
	}
	else 
		return 0;
}

string GroupPointIdentify::rotationCode(string code)
{
	if (code[0] != 'C')
		return "";
	string axisName = codeOnly(code);
	stringstream axisTypeLine;
	axisTypeLine << axisName[0] << axisName[1] << axisName[2];
	for (size_t i = 3; i < axisName.size(); i++)
	{
		if ((axisName[i] == '-')|| (axisName[i] == ' '))
		{
			break;
		}
		axisTypeLine << axisName[i];
	}
	return axisTypeLine.str();
}

int GroupPointIdentify::sImproperRotationNumber(string code)
{
	if (code[0] = 'S')
	{
		stringstream convert;
		for(size_t i = 1; i < code.size(); i++)
		{
			if (code[i] == '-')
				break;
			convert << code[i];
		}
		int nRot;
		convert >> nRot;
		return nRot;
	}
	return -1;
}

bool GroupPointIdentify::findCodeOnParentheses(
	std::string codeParentheses,
	std::string targetCode)
{
	targetCode = codeOnly(targetCode);
	stringstream convert;
	convert << codeParentheses;
	while (!convert.eof())
	{
		string tryCode;
		convert >> tryCode;
		if (tryCode == targetCode)
			return true;
	}
	return false;
}

