#include "printToFileHandler.h"


void outputHandler::initializeOutputfile(std::string outputFileName)
{
	outputStreamHandler.open(outputFileName.c_str());
}

void outputHandler::initializeOutputfile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC)
{
	std::string localName = stringA + patch::to_string(intA) + stringB + patch::to_string(intB) + stringC+".txt";
	outputStreamHandler.open(localName.c_str());
}

void outputHandler::printLineToOutputFile(std::string line)
{
	outputStreamHandler << line.c_str() << std::endl;

}
void outputHandler::printLineToOutputFile(std::string stringA, int intA)
{

	outputStreamHandler << stringA.c_str() << intA<<std::endl;

}
void outputHandler::printLineToOutputFile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC, int intC)
{

	outputStreamHandler << stringA.c_str() << intA << stringB.c_str() << intB << stringC.c_str() << intC << std::endl;


}

void outputHandler::printLineToOutputFile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC)
{


	outputStreamHandler << stringA.c_str() << intA << stringB.c_str() << intB << stringC.c_str() << std::endl;


}

void outputHandler::printLineToOutputFile(std::string stringA, int intA, std::string stringB, int intB, std::string stringC, double floatA, std::string stringD, double floatB, std::string stringE, double floatC)
{


	outputStreamHandler << stringA.c_str() << intA << stringB.c_str() << intB << stringC.c_str() << floatA << stringD.c_str() << floatB << stringE.c_str() << floatC << std::endl;


}

void outputHandler::closeOpenFile()
{
	outputStreamHandler.close();


}


