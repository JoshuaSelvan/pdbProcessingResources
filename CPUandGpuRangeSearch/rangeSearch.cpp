#include "rangeSearch.h"




void initiateRangeSearch(rangeSearchSettings RangeSearchSettings)
{
	//clock_t start, stop, start2,stop2;
	//auto start,stop;
	int totalNumberOfFiles = checkNumberOfProteinFilesInList(RangeSearchSettings.inputListFileLocation);
	AtomToNumHashTable atomReferenceTable;
	ProteinDataHandler heldProteinSets(RangeSearchSettings);
	
		

	 typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::duration<float> fsec;
    auto start =Time::now();
   auto stop = Time::now();
    auto t0 = Time::now();
    auto t1 = Time::now();
    //fsec fs = t1 - t0;
   // ms d = std::chrono::duration_cast<ms>(fs);
    //std::cout << fs.count() << "s\n";
   // std::cout << d.count() << "ms\n";




	


	t0= Time::now();
	 //start=clock();
	//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	heldProteinSets.loadAllProteinsToArrays(RangeSearchSettings.inputListFileLocation, atomReferenceTable, RangeSearchSettings);
	//std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();


	//stop=clock();
	t1 = Time::now();

	  //      print_elapsed(start2, stop2, /*"CPU KDTREE CONSTRUCTION TIME: "*/"load from esbtl time");
		print_chrono_elapsed(t0,t1, "load from esbtl time CHRONO VERSION");

	
	//  auto start2 = Time::now();
	//t0 = Time::now();
	//start2=clock();
	heldProteinSets.populateXYZarrayForGPUProcessing();
	//start2=clock();
	//
	t0 =Time::now();
	constructKdTreesOnLoadedDataOnCPU(heldProteinSets);//should this be included in the run time of searching gpu kd trees?
	//constructKdTreesOnLoadedDataOnCPUWithGPUSorting(heldProteinSets);//GPU preprocess
//constructKdTreesOnLoadedDataOnGPU(heldProteinSets);
	//GPUTreeWithCPULoad(heldProteinSets);
//	stop2 = clock();
	  //auto stop2 = Time::now();
	
//	getCountsOfAllAtomAAndBInSet(RangeSearchSettings, heldProteinSets, atomReferenceTable);
	
//	return;
	t1=Time::now();
	print_chrono_elapsed(t0, t1, "KDTREE CONSTRUCTION TIME: ");
	std::cout << std::endl;


		
//	return;
	//  auto start3 = Time::now();
	t0=Time::now();
	//clock_t start3, stop3;
	//start3 = clock();


	heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?

// Test values
//std::cout <<"OG N stats:"<<std::endl;	
//RangeSearchSettings.AtomTypeOne = "OG";
//	RangeSearchSettings.AtomTypeTwo = "N";	
//getCountsOfAllAtomAAndBInSet(RangeSearchSettings, heldProteinSets, atomReferenceTable);

//std::cout <<"CB CG stats:"<<std::endl;

//	RangeSearchSettings.AtomTypeOne = "CB";
 //       RangeSearchSettings.AtomTypeTwo = "CG";
//getCountsOfAllAtomAAndBInSet(RangeSearchSettings, heldProteinSets, atomReferenceTable);

//std::cout <<"C1' C2' stats:"<<std::endl;

//        RangeSearchSettings.AtomTypeOne = "C1'";
//        RangeSearchSettings.AtomTypeTwo = "C2'";
//getCountsOfAllAtomAAndBInSet(RangeSearchSettings, heldProteinSets, atomReferenceTable);


//return;
 // auto stop3 = Time::now();
	t1=Time::now();
//stop3 = clock();
	print_chrono_elapsed(t0,t1, "CPU SECONDARY SEARCH STRUCTURE CONSTRUCTION TIME: ");
	std::cout<<std::endl;

//	generateStatistics(RangeSearchSettings, heldProteinSets, atomReferenceTable);
	
	 if(RangeSearchSettings.searchType==-1){
                start = Time::now(); //clock();
                hybridGpuCpuSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
                stop = Time::now();//clock();
                std::cout << "RUN time for Hybrid BRUTE FORCE: ";
                print_chrono_elapsed(start, stop, "");
        }
	else if(RangeSearchSettings.searchType==0){
		start =Time::now();// clock();
		pureBruteForce(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = Time::now();//clock();
		std::cout << "RUN time for PURE CPU BRUTE FORCE: ";
		print_chrono_elapsed(start, stop, "");
	}	
	else if (RangeSearchSettings.searchType== 1){
		start = Time::now();//clock();
		//heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?
		cpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = Time::now();//clock();
		print_chrono_elapsed(start, stop, "Run time for CPU BRUTE FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 2){
		start =Time::now();// clock();
		gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = Time::now();//clock();
		print_chrono_elapsed(start, stop, "Run time for SIMPLE GPU BRUTE FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 3){
		start =Time::now();// clock();
		bruteForceSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop =Time::now();// clock();
		print_chrono_elapsed(start, stop, "Run time for SINGLE LOAD GPU BRUTE FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 4){
		start =Time::now();// clock();
		bruteForceSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = Time::now();//clock();
		print_chrono_elapsed(start, stop, "Run time for BATCH GPU Brute FORCE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 5){
		start = Time::now();//clock();
		//heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?
		cpuKdTreeRangeSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = Time::now();//clock();
		print_chrono_elapsed(start, stop, "Run time for CPU KD TREE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 6){
	
		start =Time::now();// clock();
		gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop =Time::now();// clock();
		print_chrono_elapsed(start, stop, "Run time for SIMPLE GPU KD TREE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 7){

		start = Time::now();//clock();
		kdTreeSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop =Time::now();// clock();
		print_chrono_elapsed(start, stop, "Run time for SINGLE LOAD GPU KD TREE RANGE SEARCH: ");
	} else if (RangeSearchSettings.searchType == 8)
	{

		start = Time::now();//clock();
		kdTreeSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
		stop = Time::now();//clock();
		print_chrono_elapsed(start, stop, "Run time for GPU BATCH LOAD KD TREE RANGE SEARCH: ");
	}
	else if (RangeSearchSettings.searchType == 9)
	{
		multipleRunTimer(RangeSearchSettings, heldProteinSets, atomReferenceTable);
	}
	else if (RangeSearchSettings.searchType == 10)
	{

		std::cout << "Esoteric not yet included" << std::endl;
	}

	
}



void multipleRunTimer(rangeSearchSettings RangeSearchSettings, ProteinDataHandler &heldProteinSets, AtomToNumHashTable &atomReferenceTable)
{
      typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::duration<float> fsec;

	//bruteForceSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
	//return;
	std::vector<multiRunDetailsSet> multiQuerySet;
	int numberOfSets;
	//numberOfSets=loadMultiRunDetails(multiQuerySet, "multiRunTimer.txt");
	// numberOfSets=loadMultiRunDetails(multiQuerySet, "multiRunTimerTwoAtomPairsSixMethods.txt");
//	numberOfSets=loadMultiRunDetails(multiQuerySet,"multiRunTimerDifferentPDBSizes5000.txt");
	 numberOfSets=loadMultiRunDetails(multiQuerySet,"multiRunTimer.txt");

std::cout<<"working here";
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
				auto start = Time::now();//
				auto stop = Time::now();//
				//clock_t start, stop;
				//start = clock();
				//if (multiQuerySet[i].currentRange > -1)
				//{
				//for (int y = 0; y < 5; y++)
				//{
				//	originalRangeHolder[y] = heldProteinSets.ProteinDataHolder[y].heldEntries;
				//	heldProteinSets.ProteinDataHolder[y].heldEntries = 0;
				//}
				//heldProteinSets.ProteinDataHolder[multiQuerySet[i].currentRange].heldEntries = multiQuerySet[i].entriesToSearchEachRound[p];
				//}
				 if (multiQuerySet[i].searchType == -1)
                                {
                                        start = Time::now();//clock();
                                        hybridGpuCpuSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
                                        stop = Time::now();//clock();
                                        std::cout << "RUN time for Hybrid Brute FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
                                        print_chrono_elapsed(start, stop, "");
                                }

				else if (multiQuerySet[i].searchType == 0)
				{
					start =Time::now();// clock();
					pureBruteForce(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = Time::now();//clock();
					std::cout << "RUN time for PURE CPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");
				}
				else if (multiQuerySet[i].searchType == 1)
				{
					 start =Time::now();// clock();	
					//heldProteinSets.setSecondarySearchStructure();
					cpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop =Time::now();// clock();
					std::cout << "RUN time for CPU Secondary Structure BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");
					

				
				}
				else if (multiQuerySet[i].searchType == 2)
				{
					 start =Time::now();// clock();
					gpuBruteForceRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = Time::now();//clock();
					std::cout << "RUN time for SIMPLE GPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");

					
	
				}
				else if (multiQuerySet[i].searchType == 3)
				{
					 start = Time::now();//clock();
					bruteForceSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = Time::now();//clock();
					std::cout << "Run time for SINGLE LOAD GPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");

					
				}
				else if (multiQuerySet[i].searchType == 4)
				{
					 start =Time::now();// clock();
					bruteForceSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop =Time::now();// clock();
					std::cout << "RUN TIME for BATCH LOAD GPU BRUTE FORCE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");
					
				}
				else if (multiQuerySet[i].searchType == 5)
				{
					 start =Time::now();// clock();
				//	heldProteinSets.setSecondarySearchStructure();//Should this be included in the run time of searching cpu based range searches?
					cpuKdTreeRangeSearch(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop =Time::now();// clock();
					std::cout << "RUN time for CPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");

				}
				else if (multiQuerySet[i].searchType == 6)
				{
					 start =Time::now();// clock();
					gpuKdTreeUnoptimisedRangeSearchAllLoadedSets(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = Time::now();//clock();
					std::cout << "RUN time for SIMPLE GPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");
				}
				else if (multiQuerySet[i].searchType == 7)
				{
					 start =Time::now();// clock();
					kdTreeSearchAllAvailableFilesSingleGpuLoadedEntryAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = Time::now();//clock();
					std::cout << "RUN time for SINGLE LOAD GPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");
				}
				else if (multiQuerySet[i].searchType == 8)
				{
					 start =Time::now();// clock();
					kdTreeSearchAllAvailableFilesSeveralGpuLoadedEntrysAtATime(RangeSearchSettings, heldProteinSets, atomReferenceTable);
					stop = Time::now();//clock();
					std::cout << "RUN time for BATCH LOAD GPU KDTREE: Required Proximity: " << multiQuerySet[i].proximities[j] << " Entry count: " << multiQuerySet[i].entriesToSearchEachRound[p] << " Atom pair: " << RangeSearchSettings.AtomTypeOne << "--" << RangeSearchSettings.AtomTypeTwo << std::endl;
					print_chrono_elapsed(start, stop, "");
				}

				if (multiQuerySet[i].currentRange > -1)
				{
					//for (int y = 0; y < 5; y++)
					//{
					//	heldProteinSets.ProteinDataHolder[y].heldEntries = originalRangeHolder[y];
					//}
					
				}

				
				std::cout << std::endl;
			}

			std::cout << std::endl << std::endl;
		}

	}


	return;

}

