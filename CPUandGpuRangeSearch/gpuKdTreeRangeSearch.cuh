#ifndef GPUKDTREERANGESEARCH
#define GPUKDTREERANGESEARCH
#include "cuda_runtime.h"

#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "AtomNameToNumHashTable.h"
#include <string>
#include <vector>


void gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable);
void rangeSearchCpu(atomCoords targetAtom, int elementTwo, int maxDistance, std::vector<int> &resultsVector, short *names, int *xValueArray, int *yValueArray, int*zValueArray, int *kdArray, int kdArraySize, int TreePos, int treeLevel, int* runcount, int valueArraysStartPositionOffset, int kdTreeStartPositionOffset);


//int squaredDistance3D(atomCoords targetAtom, atomCoords currentAtom);
//bool viableDirection(atomCoords targetAtom, atomCoords currentAtom, atomCoords childAtom, int treeLevel, int requiredDistance);

void calculateInitialBlocksAndThreads(int &blocks, int&threads, int maxEntrySize);
__global__ void device_side__locateElement(short * d_names, int entryNo, int soughtAtomNum, int* d_soughtAtomPositionList, int* d_soughtAtomCount, int*d_kdSearchCurrentDimensionList, int*d_kdSearchCurrentTreePos, int lengthOfChain);
__global__ void device_side_ProcessCurrentTreePositionsV3(int* d_xyzValues, short* d_Names, int*  d_kdSetArray, int* d_ViableElementAPairs, int* d_ViableElementBPairs, int* d_ViableElementPairCount, int* d_currentSearchLocations, int*d_currentSearchAElementPos, int* d_currentSearchDimensions, int* CurrentSearchCount, int SizeOfDimensionArray, int maxDistanceSquared, int elementB, int*d_nextSearchCount, int entryNum, int sizeOfKdTree, int maxDistance, int entryBeingWorkedOn);
__global__ void SetNextCountAsCurrentAndCheckFlag(int* d_currentSearchCount, int* d_nextSearchCount, int* d_completionFlag);
#endif