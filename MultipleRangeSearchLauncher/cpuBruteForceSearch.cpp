#include "cpuBruteForceSearch.h"

void cpuBruteForceRangeSearchAllLoadedSets(rangeSearchSettings &settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	int outputType = settings.resultsPrintFormat;
	outputHandler filePrinter;
	std::string printType;

	if (outputType == 3)
	{
		printType="_Summary";
	}
	else if (outputType == 4)
	{
		printType = "_Detailed";
	}

	std::cout<<std::endl << "PERFORMING CPU BRUTE FORCE RANGE SEARCH" << std::endl;

	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeOne);
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeTwo);
	int maxDistanceSquared = settings.requiredProximity*settings.requiredProximity;


	rangeSearchArrays rangeSearch;

	int allocatedSearchSize = 16384 * sizeof(int) * 4;

	rangeSearch.atomAPositionList = (int*)malloc(sizeof(int) * 16384);
	rangeSearch.atomBPositionList = (int*)malloc(sizeof(int) * 16384);
	rangeSearch.atomACount = (int*)malloc(sizeof(int) * 1);
	rangeSearch.nextSearchCount = (int*)malloc(sizeof(int) * 1);
	rangeSearch.completionFlag = (int*)malloc(sizeof(int) * 1);
	rangeSearch.MatchesCount = (int*)malloc(sizeof(int) * 1);
	rangeSearch.atomAMatches = (int*)malloc(sizeof(int) * 16384);
	rangeSearch.atomBMatches = (int*)malloc(sizeof(int) * 16384);

	int* atomBCount = (int*)malloc(sizeof(int) * 1);
	
	
	
	for (int i = 0; i < 5; i++)
	{
		if (outputType == 3 || outputType == 4)//will move into loop shortly
		{
			filePrinter.initializeOutputfile("Range_", settings.maxProteinLengthPerRange[i], "_Files_", heldProteinSets.ProteinDataHolder[i].heldEntries, printType);
		}

		if (heldProteinSets.ProteinDataHolder[i].heldEntries>0)
		{

			std::cout << std::endl;
			std::cout << "Processing Range set: " << i << std::endl;
			std::cout << "Number of present entries is: " << heldProteinSets.ProteinDataHolder[i].heldEntries << std::endl;
			if (outputType == 3 || outputType == 4)
			{
				filePrinter.printLineToOutputFile("Processing Range set: ", i);
				filePrinter.printLineToOutputFile("Number of present entries is : ", heldProteinSets.ProteinDataHolder[i].heldEntries);
				filePrinter.printLineToOutputFile("");
			}



			int numberOfProcessedFiles = 0;
			clock_t n,m;
			n=clock();


			//For every stored data set
			for (int currentEntry = 0; currentEntry < heldProteinSets.ProteinDataHolder[i].heldEntries; currentEntry++)
			{
				searchEntryInSecondaryPositionStructureForAtom(soughtAtomANumber, currentEntry, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, heldProteinSets.ProteinDataHolder[i].compositionCountsList, heldProteinSets.ProteinDataHolder[i].compositionLists, heldProteinSets.ProteinDataHolder[i].compositionPointers, rangeSearch.atomAPositionList, rangeSearch.atomACount);
				searchEntryInSecondaryPositionStructureForAtom(soughtAtomBNumber, currentEntry, heldProteinSets.ProteinDataHolder[i].MaxEntrySize, heldProteinSets.ProteinDataHolder[i].compositionCountsList, heldProteinSets.ProteinDataHolder[i].compositionLists, heldProteinSets.ProteinDataHolder[i].compositionPointers, rangeSearch.atomBPositionList, atomBCount);
				int p = 4;



				if ((rangeSearch.atomACount[0] > 0) && (atomBCount[0] > 0))
				{
					int resultsSize = 0;
					for (int p = 0; p < rangeSearch.atomACount[0]; p++)
					{

						std::vector<int> resultsVector;

						atomCoords atomACoords;
						setAtomCoords(atomACoords, heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]], heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]], heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]]);

						bruteForceSearchCpu(atomACoords, rangeSearch.atomBPositionList, atomBCount, maxDistanceSquared, heldProteinSets,resultsVector,i);

						
						for (int j = 0; j < resultsVector.size(); j++) //All atomB will be related to the atomA being searched for in the greater loop
						{
							resultsSize = resultsSize + 1;
							if (outputType == 2)
							{
								std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[rangeSearch.atomAPositionList[p]] << "\t - Pos : " << (rangeSearch.atomAPositionList[p] - heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry) << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << std::endl;
								std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry +*/ resultsVector[j]] << "\t - Pos : " << (resultsVector[j] - currentEntry*heldProteinSets.ProteinDataHolder[i].MaxEntrySize) << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry +*/ resultsVector[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + */resultsVector[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry +*/ resultsVector[j]])) / 1000) << std::endl << std::endl;
							}
							if (outputType == 4)
							{
								filePrinter.printLineToOutputFile("Atom A: " , heldProteinSets.ProteinDataHolder[i].namesSets[rangeSearch.atomAPositionList[p]] , "\t - Pos : " , (rangeSearch.atomAPositionList[p] - heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry) , "\t X: " , ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) , "\t Y: " , ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) , "\t Z: " , ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000));
								filePrinter.printLineToOutputFile("Atom B: " , heldProteinSets.ProteinDataHolder[i].namesSets[resultsVector[j]] , "\t - Pos : " , (resultsVector[j] - currentEntry*heldProteinSets.ProteinDataHolder[i].MaxEntrySize) , "\t X: " , ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[resultsVector[j]])) / 1000) , "\t Y: " , ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[resultsVector[j]])) / 1000) , "\t Z: " , ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[resultsVector[j]])) / 1000) );
								filePrinter.printLineToOutputFile("");
							}
						}


						
					}
					if (outputType == 1 || outputType == 2)
						std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: " << resultsSize << std::endl;
					else if (outputType == 3 || outputType == 4)
						filePrinter.printLineToOutputFile("Number of matches in file ", numberOfProcessedFiles, " in set ", i, "  is: ", resultsSize);
					numberOfProcessedFiles++;

				}
				else
				{
					if (outputType == 1 || outputType == 2)
						std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: 0"  << std::endl;
					else if (outputType == 3 || outputType == 4)
						filePrinter.printLineToOutputFile("Number of matches in file ", numberOfProcessedFiles, " in set ", i, "  is: 0");
					numberOfProcessedFiles++;
				}


			}
			
			m=clock();
			if (settings.debugLevel>0)
				print_elapsed(n, m, "time to process file subset: ");
		}
		filePrinter.closeOpenFile();
	}
	
}


void bruteForceSearchCpu(atomCoords atomACoords, int * atomBPositionList, int * atomBCount, int maxDistanceSquared, ProteinDataHandler heldProteinSets, std::vector<int> &resultsVector, int currentSet)
{
	for (int i = 0; i < atomBCount[0]; i++)
	{
		atomCoords atomBCoords;
		setAtomCoords(atomBCoords, heldProteinSets.ProteinDataHolder[currentSet].xCoordsSets[atomBPositionList[i]], heldProteinSets.ProteinDataHolder[currentSet].yCoordsSets[atomBPositionList[i]], heldProteinSets.ProteinDataHolder[currentSet].zCoordsSets[atomBPositionList[i]]);

		if (pow((atomACoords.xyz[0] - atomBCoords.xyz[0]), 2) + pow((atomACoords.xyz[1] - atomBCoords.xyz[1]), 2) + pow((atomACoords.xyz[2] - atomBCoords.xyz[2]), 2) <= maxDistanceSquared)
		{
			resultsVector.push_back(atomBPositionList[i]);
		}

	}
};



