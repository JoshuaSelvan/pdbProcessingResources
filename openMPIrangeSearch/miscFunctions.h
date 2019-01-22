#ifndef MISCFUNCTIONS
#define MISCFUNCTIONS


#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <vector>
#include "dataStructures.h"
//#include "RangeSearchSettings.h"
#include <chrono>

void loadRunDetails(std::string sourceFile, rangeSearchSettings &currentSettings);
void print_elapsed(clock_t start, clock_t end, std::string captionBeforeTimeDifference);
int loadMultiRunDetails(std::vector<multiRunDetailsSet> &multiSet, std::string fileName);
void resetMultiRunDetails();

void print_chrono_elapsed(std::chrono::time_point<std::chrono::high_resolution_clock>start, std::chrono::time_point<std::chrono::high_resolution_clock> end, std::string caption);

void resetInputTextFile();

void reformatFileList(std::string unformattedFile, std::string newPrefix, std::string newSuffix);

int largestSizeInArrayPos(int *theArray, int arraySize);

void switchLoadingAndProcessingSets(int &loading, int &processing);

void loadMultiBatchRangeSizes(std::string sourceFile, rangeSearchSettings &currentSettings);

void resetMultiBatchRangeSizes();


void print3ArraysToFile(int* arrayA, int* arrayB, int* arrayC, int length);

void print1ArraysToFile(int* arrayA, int length);

#endif
