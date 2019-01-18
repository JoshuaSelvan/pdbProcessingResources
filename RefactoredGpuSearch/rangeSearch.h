#ifndef CPURANGESEARCH
#define CPURANGESEARCH
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "SearchStructureConstruction.h"
#include "AtomNameToNumHashTable.h"
#include "cpuBruteForceSearch.h"
#include "gpuKdTreeRangeSearch.cuh"
#include "gpuBruteForceSearch.cuh"
#include <string>
#include <vector>



void initiateRangeSearch(rangeSearchSettings RangeSearchSettings);

void setSecondarySearchStructure(ProteinDataHandler &ProteinData);



#endif