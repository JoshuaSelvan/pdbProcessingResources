#ifndef MISCFUNCTIONS
#define MISCFUNCTIONS


#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <vector>
#include "dataStructures.h"
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
//#include "RangeSearchSettings.h"

void loadRunDetails(std::string sourceFile, rangeSearchSettings &currentSettings);
void print_elapsed(clock_t start, clock_t end, std::string captionBeforeTimeDifference);

void print_chrono_elapsed(std::chrono::time_point<std::chrono::high_resolution_clock>start, std::chrono::time_point<std::chrono::high_resolution_clock> end, std::string caption);

int loadMultiRunDetails(std::vector<multiRunDetailsSet> &multiSet, std::string fileName);
void resetMultiRunDetails();
//void setIntArrayValues(int *array, int value, int lengthOfArray);
//void setShortArrayValues(short *array, int value, int lengthOfArray);

void resetInputTextFile();

void reformatFileList(std::string unformattedFile, std::string newPrefix, std::string newSuffix, std::string newFileName);

int largestSizeInArrayPos(int *theArray, int arraySize);

void switchLoadingAndProcessingSets(int &loading, int &processing);

void loadMultiBatchRangeSizes(std::string sourceFile, rangeSearchSettings &currentSettings);

void resetMultiBatchRangeSizes();

#endif
