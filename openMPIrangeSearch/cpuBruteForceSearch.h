#ifndef CPUBRUTEFORCERANGESEARCH
#define CPUBRUTEFORCERANGESEARCH
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "AtomNameToNumHashTable.h"
#include "SearchStructureConstruction.h"
#include "dataStructures.h"
#include <string>
#include <vector>
#include <cmath>
#include "printToFileHandler.h"

void cpuBruteForceRangeSearchAllLoadedSets( rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);


void bruteForceSearchCpu(atomCoords atomACoords, int * atomBPositionList, int * atomBCount, int maxDistanceSquared, ProteinDataHandler heldProteinSets, std::vector<int>& resultsVector, int currentSet);


void pureBruteForce(rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);

#endif
