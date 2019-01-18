#ifndef PROTEINLOADER
#define PROTEINLOADER

#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
#include "AtomNameToNumHashTable.h"
#include "miscFunctions.h"
//#include "SearchStructureConstruction.h"
#include "dataStructures.h"
#include "ESBTL_ProteinHandler.h"
#include <string>
#include <assert.h>

#include <algorithm>    // std::fill_n



int checkNumberOfProteinFilesInList(std::string inputFileLocationList); //checks how many files are listed in an inputFile;

class ProteinDataHandler{

public:
    
	multipleProteinSet ProteinDataHolder[5];
	ProteinDataHandler();
	
	void  createSubListsContainingValidFilesForEachRangeOnly(rangeSearchSettings &currentSettings);
	int determineRangeOfSingleProteinFile(std::string proteinFileLocation, rangeSearchSettings &currentSettings);


    ProteinDataHandler(rangeSearchSettings currentSettings); 
	
	void loadAllProteinsToArrays(std::string inputFileLocationsList, AtomToNumHashTable &atomReferenceTable, rangeSearchSettings &currentSettings); //returns number of atoms found in current protein
	
	
	int loadSingleProteinDetailsIntoArrays(std::string proteinFileLocation, AtomToNumHashTable &atomReferenceTable, rangeSearchSettings &currentSettings);//Places a single set of protein data into an availalbe space in the multipleProteinSet.
	void populateXYZarrayForGPUProcessing();


	void setSecondarySearchStructure();


};


#endif
