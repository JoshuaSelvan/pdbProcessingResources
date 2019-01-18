#ifndef MISCFUNCTIONS
#define MISCFUNCTIONS


#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <vector>
#include "dataStructures.h"
//#include "RangeSearchSettings.h"

void loadRunDetails(std::string sourceFile, rangeSearchSettings &currentSettings);
void print_elapsed(clock_t start, clock_t end, std::string captionBeforeTimeDifference);
int loadMultiRunDetails(std::vector<multiRunDetailsSet> &multiSet, std::string fileName);
void resetMultiRunDetails();
//void setIntArrayValues(int *array, int value, int lengthOfArray);
//void setShortArrayValues(short *array, int value, int lengthOfArray);

void resetInputTextFile();

void reformatFileList(std::string unformattedFile, std::string newPrefix, std::string newSuffix);

#endif
