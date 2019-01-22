#include "rangeSearchHandler.h"




void initiateRangeSearch(rangeSearchSettings RangeSearchSettings,int personalFilesToProcess, int personalStartingPointInList)
{
	typedef std::chrono::high_resolution_clock Time;
    	typedef std::chrono::milliseconds ms;
    	typedef std::chrono::duration<float> fsec;

	auto start = Time::now();
        auto stop = Time::now();


	AtomToNumHashTable atomReferenceTable;
	ProteinDataHandler heldProteinSets(RangeSearchSettings);
	

	 start = Time::now();

	heldProteinSets.loadAllProteinsToArrays(RangeSearchSettings.inputListFileLocation, atomReferenceTable, RangeSearchSettings, personalFilesToProcess, personalStartingPointInList);
	
	 stop = Time::now();
 print_chrono_elapsed(start, stop, "ESBTL LOADING TIME: ");

	heldProteinSets.populateXYZarrayForGPUProcessing();
	 start = Time::now();	
//start2=clock();
	constructKdTreesOnLoadedDataOnCPU(heldProteinSets);
	//constructKdTreesOnLoadedDataOnCPUWithGPUSorting(heldProteinSets);
	//constructKdTreesOnLoadedDataOnGPU(heldProteinSets);
	//GPUTreeWithCPULoad(heldProteinSets);
	//stop2 = clock();
	stop = Time::now();
	print_chrono_elapsed(start, stop, "CPU KDTREE CONSTRUCTION TIME: ");
	//std::cout << std::endl;
	

	//clock_t start3, stop3;
	//start3 = clock();
	 start = Time::now();
	heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?
	 stop = Time::now();
	print_chrono_elapsed(start, stop, "CPU SECONDARY SEARCH STRUCTURE CONSTRUCTION TIME: ");


	//generateStatistics(RangeSearchSettings, heldProteinSets, atomReferenceTable);
	//return;

	if (RangeSearchSettings.searchType == 0){
	 start = Time::now();	
	//start = clock();
		pureBruteForce(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		//stop = clock();
		 stop = Time::now();
		std::cout << "RUN time for PURE CPU BRUTE FORCE: Required Proximity: ";
		print_chrono_elapsed(start, stop, "");
	}
	else if (RangeSearchSettings.searchType== 1){
		//start = clock();
		 start = Time::now();
		cpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		//stop = clock();
		 stop = Time::now();
		print_chrono_elapsed(start, stop, "Run time for CPU BRUTE FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 5){
		//start = clock();
		 start = Time::now();
		//heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?
		cpuKdTreeRangeSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		//stop = clock();
		 stop = Time::now();
		print_chrono_elapsed(start, stop, "Run time for CPU KD TREE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 9)
	{
		multipleRunTimer(RangeSearchSettings, heldProteinSets, atomReferenceTable);
	}
	else if (RangeSearchSettings.searchType == 10)
	{

		std::cout << "Esoteric not included" << std::endl;
	}

	
}



void multipleRunTimer(rangeSearchSettings RangeSearchSettings, ProteinDataHandler &heldProteinSets, AtomToNumHashTable &atomReferenceTable)
{
	std::vector<multiRunDetailsSet> multiQuerySet;
	int numberOfSets;
	numberOfSets=loadMultiRunDetails(multiQuerySet, "multiRunTimer.txt");
	int originalRangeHolder[5];

	for (int i = 0; i < numberOfSets; i++)
	{

		RangeSearchSettings.AtomTypeOne = multiQuerySet[i].atom1;
		RangeSearchSettings.AtomTypeTwo = multiQuerySet[i].atom2;
		for (int j = 0; j < multiQuerySet[i].proximitiesCount; j++)
		{
			RangeSearchSettings.requiredProximity = multiQuerySet[i].proximities[j]*1000;
			for (int p = 0; p < multiQuerySet[i].entriesToSearchEachRound.size(); p++)
			{
				clock_t start, stop;
				start = clock();
				if (multiQuerySet[i].currentRange > -1)
				{
					for (int y = 0; y < 5; y++)
					{
						originalRangeHolder[y] = heldProteinSets.ProteinDataHolder[y].heldEntries;
						heldProteinSets.ProteinDataHolder[y].heldEntries = 0;
					}
					heldProteinSets.ProteinDataHolder[multiQuerySet[i].currentRange].heldEntries = multiQuerySet[i].entriesToSearchEachRound[p];
				}

				if (multiQuerySet[i].searchType == 0)
				{
					pureBruteForce(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for PURE CPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");
				}
				else if (multiQuerySet[i].searchType == 1)
				{
					
					//heldProteinSets.setSecondarySearchStructure();
					cpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for CPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");
					

				
				}
				else if (multiQuerySet[i].searchType == 2)
				{
			//		gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for SIMPLE GPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");

					
	
				}
				else if (multiQuerySet[i].searchType == 3)
				{
			//		bruteForceSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "Run time for SINGLE LOAD GPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");

					
				}
				else if (multiQuerySet[i].searchType == 4)
				{

			//		bruteForceSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN TIME for BATCH LOAD GPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");
					
				}
				else if (multiQuerySet[i].searchType == 5)
				{
					//heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?
			//		cpuKdTreeRangeSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for CPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");

				}
				else if (multiQuerySet[i].searchType == 6)
				{
			//		gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for SIMPLE GPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");
				}
				else if (multiQuerySet[i].searchType == 7)
				{
			//		kdTreeSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for SINGLE LOAD GPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");
				}
				else if (multiQuerySet[i].searchType == 8)
				{
			//		kdTreeSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = clock();
					std::cout << "RUN time for BATCH LOAD GPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_elapsed(start, stop, "");
				}


				if (multiQuerySet[i].currentRange > -1)
				{
					for (int y = 0; y < 5; y++)
					{
						heldProteinSets.ProteinDataHolder[y].heldEntries = originalRangeHolder[y];
					}
					
				}

				
				std::cout << std::endl;
			}

			std::cout << std::endl << std::endl;
		}

	}


	return;

}

