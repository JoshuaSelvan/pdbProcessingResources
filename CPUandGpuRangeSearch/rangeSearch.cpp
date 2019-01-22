#include "rangeSearch.h"




void initiateRangeSearch(rangeSearchSettings RangeSearchSettings)
{
	clock_t start, stop;
	int totalNumberOfFiles = checkNumberOfProteinFilesInList(RangeSearchSettings.inputListFileLocation);
	AtomToNumHashTable atomReferenceTable;
	ProteinDataHandler heldProteinSets(RangeSearchSettings);
	
	

	
	



	
	heldProteinSets.loadAllProteinsToArrays(RangeSearchSettings.inputListFileLocation, atomReferenceTable, RangeSearchSettings);
	constructKdTreesOnLoadedDataOnCPU(heldProteinSets);
	heldProteinSets.setSecondarySearchStructure();
	//heldProteinSets.populateXYZarrayForGPUProcessing();

	
	if (RangeSearchSettings.searchType== 1){
		start = clock();
		cpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable );
		stop = clock();
		print_elapsed(start, stop, "Run time for CPU BRUTE FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 2){
		start = clock();
		cpuKdTreeRangeSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = clock();
		print_elapsed(start, stop, "Run time for CPU KD TREE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 3){
		start = clock();
		gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = clock();
		print_elapsed(start, stop, "Run time for GPU BRUTE FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 4){
		start = clock();
		gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = clock();
		print_elapsed(start, stop, "Run time for GPU KD TREE RANGE SEARCH: ");
	}else{
		multipleRunTimer(RangeSearchSettings, heldProteinSets, atomReferenceTable);



	}
	

	//gpuKdTreeRangeSearchAllLoadedSets(RangeSearchSettings.AtomTypeOne, RangeSearchSettings.AtomTypeTwo, RangeSearchSettings.requiredProximity, heldProteinSets, atomReferenceTable);
}



void multipleRunTimer(rangeSearchSettings RangeSearchSettings, ProteinDataHandler &heldProteinSets, AtomToNumHashTable &atomReferenceTable)
{
	std::vector<multiRunDetailsSet> multiQuerySet;
	int numberOfSets;
	numberOfSets=loadMultiRunDetails(multiQuerySet, "multiRunTimer.txt");

	for (int i = 0; i < numberOfSets; i++)
	{

		RangeSearchSettings.AtomTypeOne = multiQuerySet[i].atom1;
		RangeSearchSettings.AtomTypeTwo = multiQuerySet[i].atom2;
		for (int j = 0; j < multiQuerySet[i].proximities.size(); j++)
		{
			RangeSearchSettings.requiredProximity = multiQuerySet[i].proximities[j]*1000;
			for (int p = 0; p < multiQuerySet[i].entriesToSearchEachRound.size(); p++)
			{
				clock_t start, stop;
				start = clock();
				heldProteinSets.ProteinDataHolder[multiQuerySet[i].currentRange].heldEntries = multiQuerySet[i].entriesToSearchEachRound[p];

				if (multiQuerySet[i].searchType == 0)
				{
					cpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN TIME for CPU BRUTE FORCE: "<<std::endl<<"Required Proximity: " << multiQuerySet[i].proximities[j] <<std::endl<< "Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] <<std::endl<<  "Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop,"Full run time: ");
				}
				else if (multiQuerySet[i].searchType == 1)
				{
					cpuKdTreeRangeSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN TIME for CPU KdTree: "<<std::endl<<"Required Proximity: " << multiQuerySet[i].proximities[j] <<std::endl<< "Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] <<std::endl<<  "Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "Full run time: ");
				}
				else if (multiQuerySet[i].searchType == 2)
				{
					gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN TIME for GPU BRUTE FORCE: "<<std::endl<<"Required Proximity: " << multiQuerySet[i].proximities[j] <<std::endl<<  "Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] <<std::endl<<  "Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "Full run time: ");
				}
				else if (multiQuerySet[i].searchType == 3)
				{
					gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN TIME for GPU KdTree FORCE: "<<std::endl<<"Required Proximity: " << multiQuerySet[i].proximities[j] <<std::endl<<  "Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] <<std::endl<<  "Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "Full run time: ");
				}


				
				std::cout << std::endl;
			}

			std::cout << std::endl << std::endl;
		}

	}


	return;

}



/*
float proximities[4] = { 0.5, 1, 4, 8 };
int entriesToSearch[9] = { 1000, 2000, 5000, 7000 };
int currentRange = 0;


RangeSearchSettings.AtomTypeOne = "C1'";
RangeSearchSettings.AtomTypeTwo = "C2'";
for (int i = 0; i < 4; i++)
{
	RangeSearchSettings.requiredProximity = proximities[i];

	for (int j = 0; j < 4; j++)
	{
		heldProteinSets.ProteinDataHolder[currentRange].heldEntries = entriesToSearch[j];
		start = clock();
		gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = clock();
		std::cout << "RUN TIME for GPU BRUTE FORCE: Required Proximity: " << proximities[i] << " Entry count: " << entriesToSearch[j] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
		print_elapsed(start, stop, "");
		std::cout << std::endl;

	}


	std::cout << std::endl << std::endl;
}

RangeSearchSettings.AtomTypeOne = "CB";
RangeSearchSettings.AtomTypeTwo = "CG";
for (int i = 0; i < 4; i++)
{
	RangeSearchSettings.requiredProximity = proximities[i];

	for (int j = 0; j < 4; j++)
	{
		heldProteinSets.ProteinDataHolder[currentRange].heldEntries = entriesToSearch[j];
		start = clock();
		gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = clock();
		std::cout << "RUN TIME for GPU BRUTE FORCE: Required Proximity: " << proximities[i] << " Entry count: " << entriesToSearch[j] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
		print_elapsed(start, stop, "");
		std::cout << std::endl;

	}



}*/
