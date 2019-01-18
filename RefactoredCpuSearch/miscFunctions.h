#ifndef MISCFUNCTIONS
#define MISCFUNCTIONS


#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include "dataStructures.h"
//#include "RangeSearchSettings.h"

void loadRunDetails(std::string sourceFile, rangeSearchSettings &currentSettings);
void print_elapsed(clock_t start, clock_t end, std::string captionBeforeTimeDifference);

void setIntArrayValues(int *array, int value, int lengthOfArray);
void setShortArrayValues(short *array, int value, int lengthOfArray);

void resetInputTextFile();

#endif
