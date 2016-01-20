#ifndef MYEXCEPTIONS_H
#define MYEXCEPTIONS_H

#include <exception>
#include <string>

class MyExceptions : public std::exception
{
public:
	MyExceptions(int iError_in);
	~MyExceptions();

	const char * what() const throw ();

private:
	int iError;


};

#endif
