#ifndef CPUBRUTEFORCERANGESEARCH
#define CPUBRUTEFORCERANGESEARCH
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "AtomNameToNumHashTable.h"
#include "SearchStructureConstruction.h"
#include <string>
#include <vector>
#include <cmath>

void bruteForceRangeSearchAllLoadedSets(std::string AtomA, std::string atomB, int requiredProximity, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);
//void bruteForceSearchCpu(atomCoords targetAtom, int elementTwo, int maxDistance, std::vector<int> &resultsVector, short *names, int *xValueArray, int *yValueArray, int*zValueArray, int *kdArray, int kdArraySize, int TreePos, int treeLevel, int* runcount, int valueArraysStartPositionOffset, int kdTreeStartPositionOffset);


void bruteForceSearchCpu(atomCoords atomACoords, int * atomBPositionList, int * atomBCount, int maxDistanceSquared, ProteinDataHandler heldProteinSets, std::vector<int>& resultsVector, int currentSet);

#endif
