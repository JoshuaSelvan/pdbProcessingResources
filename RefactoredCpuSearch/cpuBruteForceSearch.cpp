#include "cpuBruteForceSearch.h"

void bruteForceRangeSearchAllLoadedSets(std::string atomA, std::string atomB, int requiredProximity, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	std::cout << "PERFORMING BRUTE FORCE SEARCH" << std::endl;

	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(atomA);//heldProteinSets.ProteinDataHolder[0].namesSets[0];//
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(atomB);//heldProteinSets.ProteinDataHolder[0].namesSets[1];//
	int maxDistanceSquared = requiredProximity*requiredProximity;


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
		if (heldProteinSets.ProteinDataHolder[i].heldEntries>0)
		{
			std::cout << std::endl;
			std::cout << "Processing Range set: " << i << std::endl;
			std::cout << "Number of present entries is: " << heldProteinSets.ProteinDataHolder[i].heldEntries << std::endl;

			int numberOfProcessedFiles = 0;

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
							/*outputFile*/ //std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[rangeSearch.atomAPositionList[p]] << "\t - Pos : " << (rangeSearch.atomAPositionList[p] - heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry) << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << std::endl;
							/*outputFile*/ //std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry +*/ resultsVector[j]] << "\t - Pos : " << (resultsVector[j]-currentEntry*heldProteinSets.ProteinDataHolder[i].MaxEntrySize) << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry +*/ resultsVector[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + */resultsVector[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[/*heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry +*/ resultsVector[j]])) / 1000) << std::endl << std::endl;
						
							
						
						}


						//resultsSize = resultsSize + resultsVector.size();

					}
					std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: " << resultsSize << std::endl;
					numberOfProcessedFiles++;

				}
				else
				{
					std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: 0"  << std::endl;
					numberOfProcessedFiles++;
				}


			}
			//	outputFile.close();
			
		}
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



