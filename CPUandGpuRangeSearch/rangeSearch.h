#ifndef CPURANGESEARCH
#define CPURANGESEARCH
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "SearchStructureConstruction.h"
#include "AtomNameToNumHashTable.h"
#include "cpuBruteForceSearch.h"
#include "gpuKdTreeRangeSearch.cuh"
#include "gpuBruteForceSearch.cuh"
#include "cpuKdTreeSearch.h"
//#include "StatisticsGenerator.h"
#include <string>
#include <vector>
#include "constructKdTreesOnLoadedDataWithGPU.cuh"

#include <chrono>

//#include "chrono_io"





void initiateRangeSearch(rangeSearchSettings RangeSearchSettings);
void multipleRunTimer(rangeSearchSettings RangeSearchSettings, ProteinDataHandler &heldProteinSets, AtomToNumHashTable &atomReferenceTable);

//void setSecondarySearchStructure(ProteinDataHandler &ProteinData);



#endif
