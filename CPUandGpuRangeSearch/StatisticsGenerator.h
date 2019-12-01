#ifndef STATISTICSGENERATOR
#define STATISTICSGENERATOR
#include "ESBTL_ProteinHandler.h"
#include "dataStructures.h"
#include "AtomNameToNumHashTable.h"
#include "SearchStructureConstruction.h"
#include <string>
#include <vector>
#include <cmath>
#include "printToFileHandler.h"

void generateStatistics(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable);


#endif

