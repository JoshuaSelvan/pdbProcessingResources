#ifndef GPUBRUTEFORCERANGESEARCH
#define GPUBRUTEFORCERANGESEARCH
#include "cuda_runtime.h"

#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "AtomNameToNumHashTable.h"
#include "printToFileHandler.h"
#include "miscFunctions.h"
#include "SearchStructureConstruction.h"
#include <string>
#include <vector>


void gpuBruteForceRangeSearchAllLoadedSets(/*std::string AtomA, std::string atomB, int requiredProximity, */ rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);


void calculateInitialBlocksAndThreads(int &blocks, int&threads, int maxEntrySize);
__global__ void device_side__locateElement(short * d_names, int entryNo, int soughtAtomNum, int* d_soughtAtomPositionList, int* d_soughtAtomCount, int*d_kdSearchCurrentDimensionList, int*d_kdSearchCurrentTreePos, int lengthOfChain);



void bruteForceSearchPreLoadedArraySets(short * d_namesSet, int*d_xValsSet, int*d_yValsSet, int*d_zValsSet, int MaximumlengthOfChains, short elementA, short  elementB, int numOfEntries, rangeSearchSettings& settings, int * h_xValsSets, int * h_yValSets, int *h_zValSets, short* h_names, int currentSet);
__global__ void DeviceLoadedArrays_SingleProtein_CompleteBruteForceSearch(short * d_namesSet, int* d_elementAList, int* d_elementBList, short atomA, short atomB, int * d_aCount, int * d_bCount, int standardizedEntrySize, int entryNumber, int concurrentThreads, int*d_resultsAList, int*d_resultsBList, int*d_resultsCount, int* d_xValsSet, int* d_yValsSet, int* d_zValsSet, int requiredProximity);


__global__ void DeviceLoadedArrays_SingleProtein_LocateElements(short * d_namesSet, int* d_elementAList, int* d_elementBList, short atomA, short atomB, int * d_aCount, int * d_bCount, int standardizedEntrySize, int entryNumber, int concurrentThreads, int*d_resultsAList, int*d_resultsBList, int*d_resultsCount, int* d_xValsSet, int* d_yValsSet, int* d_zValsSet, int requiredProximity);
__global__ void DeviceLoadedArrays_SingleProtein_BruteForceSearch(short * d_namesSet, int* d_elementAList, int* d_elementBList, short atomA, short atomB, int * d_aCount, int * d_bCount, int standardizedEntrySize, int entryNumber, int concurrentThreads, int*d_resultsAList, int*d_resultsBList, int*d_resultsCount, int* d_xValsSet, int* d_yValsSet, int* d_zValsSet, int requiredProximity);


void bruteForceSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);
void bruteForceSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);
void hybridGpuCpuSearch(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);
#endif
