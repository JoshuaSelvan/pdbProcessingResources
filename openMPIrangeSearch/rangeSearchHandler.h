#ifndef CPURANGESEARCH
#define CPURANGESEARCH
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "SearchStructureConstruction.h"
#include "AtomNameToNumHashTable.h"
#include "cpuBruteForceSearch.h"

#include "cpuKdTreeSearch.h"
//#include "StatisticsGenerator.h"
#include <string>
#include <vector>

#include <chrono>


void initiateRangeSearch(rangeSearchSettings RangeSearchSettings,int personalFilesToProcess, int personalStartingPointInList);
void multipleRunTimer(rangeSearchSettings RangeSearchSettings, ProteinDataHandler &heldProteinSets, AtomToNumHashTable &atomReferenceTable);

//void setSecondarySearchStructure(ProteinDataHandler &ProteinData);



#endif
