#ifndef MYEXCEPTIONS_H
#define MYEXCEPTIONS_H

#include <exception>
#include <string>

class MyExceptions : public std::exception
{
public:
	MyExceptions(int iError_in);
	MyExceptions(std::string custom);
	~MyExceptions();

	const char * what() const throw ();

private:
	int iError;
	std::string customMessage;

};

#endif
