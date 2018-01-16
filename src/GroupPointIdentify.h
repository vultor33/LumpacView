#ifndef GROUPPOINTIDENTIFY_H
#define GROUPPOINTIDENTIFY_H

#include <string>
#include <vector>

/*format
rot     - "C2-9  ( C4-1-p C4-1-m C2-1 C3-1-p C3-1-m C3-3-p C3-3-m C2-8 )" - perpendicular in parentheses
rot     - "C3-1-p"  (only C2 have to show perpendicular)
plane   -  "P-9  ( C2-4 )"  -> rotations perpendicular in parentheses
inv     -  "Inv"
S       -  "S6-1-p"
*/
class GroupPointIdentify
{
public:
	GroupPointIdentify();

	~GroupPointIdentify();

	std::string findGroupPoint(std::vector<std::string> & allSymmetries);

	std::string codeOnly(std::string wholeCode);





private:
	int findPrincipalAxis(std::vector<std::string> & allSymmetries);

	bool findPlane(std::vector<std::string> & allSymmetries);

	bool findInversion(std::vector<std::string> & allSymmetries);

	int findOtherCnAxis(
		std::vector<std::string> & allSymmetries,
		int princAxis);

	bool findnPlanes(
		std::vector<std::string> & allSymmetries,
		int princAxis);

	bool findnC2Perpendicular(
		std::vector<std::string> & allSymmetries,
		int princAxis);

	bool findPlanePerpendicular(
		std::vector<std::string> & allSymmetries,
		int princAxis);

	bool find2nSimproper(
		std::vector<std::string> & allSymmetries,
		int princAxis);


	
	int rotationNumber(std::string code);

	int sImproperRotationNumber(std::string code);

	std::string rotationCode(std::string code);

	bool findCodeOnParentheses(
		std::string codeParentheses,
		std::string targetCode);

};




#endif