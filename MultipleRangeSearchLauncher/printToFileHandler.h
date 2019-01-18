#ifndef PRINTTOFILEHANDLER
#define PRINTTOFILEHANDLER


#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
namespace patch
{
	template < typename T > std::string to_string(const T& n)
	{
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

#include <iostream>


class outputHandler{

public:

	std::ofstream outputStreamHandler;
	void initializeOutputfile(std::string outputFileName);
	void initializeOutputfile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC);
	void printLineToOutputFile(std::string line);
	void closeOpenFile();
	void printLineToOutputFile(std::string stringA, int intA);
	void printLineToOutputFile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC, int intC);
	void printLineToOutputFile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC, double floatA, std::string stringD, double floatB, std::string stringE, double floatC);
	void printLineToOutputFile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC);
	

};
#endif
