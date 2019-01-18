#ifndef CPUKDTREERANGESEARCH
#define CPUKDTREERANGESEARCH
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "AtomNameToNumHashTable.h"
#include "SearchStructureConstruction.h"
#include <string>
#include <vector>
#include <cmath>

void rangeSearchAllLoadedSets(std::string AtomA, std::string atomB, int requiredProximity, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);
void rangeSearchCpu(atomCoords targetAtom, int elementTwo, int maxDistance, std::vector<int> &resultsVector, short *names, int *xValueArray, int *yValueArray, int*zValueArray, int *kdArray, int kdArraySize, int TreePos, int treeLevel, int* runcount, int valueArraysStartPositionOffset, int kdTreeStartPositionOffset);


int squaredDistance3D(atomCoords targetAtom, atomCoords currentAtom);
bool viableDirection(atomCoords targetAtom, atomCoords currentAtom, atomCoords childAtom, int treeLevel, int requiredDistance);

#endif
