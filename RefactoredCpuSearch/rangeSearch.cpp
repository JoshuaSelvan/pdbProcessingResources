#include "rangeSearch.h"




void initiateCpuRangeSearch(rangeSearchSettings RangeSearchSettings)
{

	


	int totalNumberOfFiles = checkNumberOfProteinFilesInList(RangeSearchSettings.inputListFileLocation);
	AtomToNumHashTable atomReferenceTable;
	ProteinDataHandler heldProteinSets(RangeSearchSettings);


	heldProteinSets.loadAllProteinsToArrays(RangeSearchSettings.inputListFileLocation, atomReferenceTable);
	constructKdTreesOnLoadedData(heldProteinSets);
	setSecondarySearchStructure(heldProteinSets);
    

	clock_t start,stop;
	start=clock();
	
	//Select range search method.
	if (RangeSearchSettings.searchType==1)
	{
		bruteForceRangeSearchAllLoadedSets(RangeSearchSettings.AtomTypeOne, RangeSearchSettings.AtomTypeTwo, RangeSearchSettings.requiredProximity, heldProteinSets, atomReferenceTable);
	}
	else if (RangeSearchSettings.searchType == 2)
	{
		rangeSearchAllLoadedSets(RangeSearchSettings.AtomTypeOne, RangeSearchSettings.AtomTypeTwo, RangeSearchSettings.requiredProximity, heldProteinSets, atomReferenceTable);
	}
	else if ((RangeSearchSettings.searchType == 3) || (RangeSearchSettings.searchType == 4)){
		std::cout << "Requested search option  not available" << std::endl;
	}
	else{
		std::cout << "Requested search option doesn't seem to exist" << std::endl;
	}

	stop=clock();


	print_elapsed(start,stop,"time to perform search");
}



void setSecondarySearchStructure(ProteinDataHandler &ProteinData)
{
	//Null all secondary Structure arrays;
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < ProteinData.ProteinDataHolder[i].heldEntries; j++)
		{
			ProteinData.ProteinDataHolder[i].compositionCountsList[j] = 0; //holds the number of atom types in each protein
		}
	}

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < ProteinData.ProteinDataHolder[i].heldEntries*ProteinData.ProteinDataHolder[i].KdTreeSize; j++)
		{
			ProteinData.ProteinDataHolder[i].compositionLists[j] = 0;//Holds pointers to the positions of all the atoms in each protein set
			ProteinData.ProteinDataHolder[i].compositionPointers[j] = 0;//Will hold the names of each aotm type present followed by the position in the composition list where those atoms start
		}
	}

	//Construct secondary structures 
	for (int i = 0; i<5; i++)
	{
		formatSecondaryPositionStructure(ProteinData.ProteinDataHolder[i].namesSets, ProteinData.ProteinDataHolder[i].MaxEntrySize, ProteinData.ProteinDataHolder[i].heldEntries, ProteinData.ProteinDataHolder[i].proteinLengthCounts, ProteinData.ProteinDataHolder[i].compositionCountsList, ProteinData.ProteinDataHolder[i].compositionLists, ProteinData.ProteinDataHolder[i].compositionPointers);
	}



}

