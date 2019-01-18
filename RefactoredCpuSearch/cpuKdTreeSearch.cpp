#include "cpuKdTreeSearch.h"

void rangeSearchAllLoadedSets(std::string atomA, std::string atomB, int requiredProximity, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{
	std::cout << "PERFORMING KD TREE RANGE SEARCH" << std::endl;
	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(atomA);//heldProteinSets.ProteinDataHolder[0].namesSets[0];
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(atomB);//heldProteinSets.ProteinDataHolder[0].namesSets[1];
	int maxDistanceSquared = requiredProximity*requiredProximity;


	rangeSearchArrays rangeSearch;

	int allocatedSearchSize = 16384 * sizeof(int) * 4;

	rangeSearch.atomAPositionList = (int*)malloc(sizeof(int) * 16384);
	rangeSearch.atomACount = (int*)malloc(sizeof(int) * 1);
	//rangeSearch.nextSearchCount = (int*)malloc(sizeof(int) * 1);
	rangeSearch.completionFlag = (int*)malloc(sizeof(int) * 1);
	//rangeSearch.MatchesCount = (int*)malloc(sizeof(int) * 1);
	rangeSearch.atomAMatches = (int*)malloc(sizeof(int) * 16384);
	rangeSearch.atomBMatches = (int*)malloc(sizeof(int) * 16384);

	rangeSearch.atomACount[0] = 0;
	for (int i = 0; i < 16384; i++)
	{
		rangeSearch.atomAMatches[i] = -1;
		rangeSearch.atomBMatches[i] = -1;
		rangeSearch.atomAPositionList[i] = -1;
	}

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
				if (rangeSearch.atomACount[0] > 0)
				{
					int resultsSize = 0;
					for (int p = 0; p < rangeSearch.atomACount[0]; p++)
					{

						std::vector<int> resultsVector;
						int *runCount = (int*)malloc(sizeof(int));
						runCount[0] = 0;

						atomCoords atomACoords;
						setAtomCoords(atomACoords, heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]], heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]], heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]]);

						rangeSearchCpu(atomACoords, soughtAtomBNumber, requiredProximity, resultsVector, heldProteinSets.ProteinDataHolder[i].namesSets, heldProteinSets.ProteinDataHolder[i].xCoordsSets, heldProteinSets.ProteinDataHolder[i].yCoordsSets, heldProteinSets.ProteinDataHolder[i].zCoordsSets, heldProteinSets.ProteinDataHolder[i].kdTrees, heldProteinSets.ProteinDataHolder[i].KdTreeSize, 0, 0, runCount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry, heldProteinSets.ProteinDataHolder[i].KdTreeSize*currentEntry);


						for (int j = 0; j < resultsVector.size(); j++) //All atomB will be related to the atomA being searched for in the greater loop
						{
							resultsSize = resultsSize + 1;
							/*outputFile*/ //std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[rangeSearch.atomAPositionList[p]] << "\t - Pos : " << (rangeSearch.atomAPositionList[p] - heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry) << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << std::endl;
							/*outputFile*/ //std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]] << "\t - Pos : " << resultsVector[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000) << std::endl << std::endl;
						}


						//resultsSize=resultsSize+resultsVector.size();

					}
					std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: " << resultsSize << std::endl;
					numberOfProcessedFiles++;

				}
				else
				{
					std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: 0" << std::endl;
					numberOfProcessedFiles++;
				}


			}
			//	outputFile.close();
			//std::cout << "finished range searches for set: " << i << std::endl;
		}
	}
}

void rangeSearchCpu(atomCoords targetAtom, int elementTwo, int maxDistance, std::vector<int> &resultsVector, short *names, int *xValueArray, int *yValueArray, int*zValueArray, int *kdArray, int kdArraySize, int TreePos, int treeLevel, int* runcount, int valueArraysStartPositionOffset, int kdTreeStartPositionOffset)
{
	atomCoords currentAtom;
	atomCoords rightChildAtom;
	atomCoords leftChildAtom;
	runcount[0] = runcount[0]++;
	setAtomCoords(currentAtom, xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]]);

	int currentElement = names[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + TreePos]];
	int leftChildNode = TreePos * 2 + 1;
	int rightChildNode = TreePos * 2 + 2;
	float squaredDistanceFromAtomAToLocation = squaredDistance3D(targetAtom, currentAtom);





	if (kdArray[kdTreeStartPositionOffset + TreePos] == -1)
		return;

	if (squaredDistanceFromAtomAToLocation < pow(maxDistance, 2)) //If the distance between the desired atom and the current atom is within the desired range
	{

		if (currentElement == (elementTwo))
			resultsVector.push_back(kdArray[kdTreeStartPositionOffset + TreePos]);

		if ((leftChildNode < kdArraySize) && (kdArray[kdTreeStartPositionOffset + leftChildNode] != -1))
			rangeSearchCpu(targetAtom, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, leftChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);

		if ((rightChildNode < kdArraySize) && (kdArray[kdTreeStartPositionOffset + rightChildNode] != -1))
			rangeSearchCpu(targetAtom, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, rightChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);
	}
	else
	{
		setAtomCoords(rightChildAtom, xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + rightChildNode]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + rightChildNode]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + rightChildNode]]);

		setAtomCoords(leftChildAtom, xValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + leftChildNode]], yValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + leftChildNode]], zValueArray[valueArraysStartPositionOffset + kdArray[kdTreeStartPositionOffset + leftChildNode]]);


		if ((leftChildNode < kdArraySize) && (kdArray[kdTreeStartPositionOffset + leftChildNode] != -1) && (viableDirection(targetAtom, currentAtom, leftChildAtom, treeLevel, maxDistance) == true))
			rangeSearchCpu(targetAtom, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, leftChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);

		if ((rightChildNode < kdArraySize) && (kdArray[kdTreeStartPositionOffset + rightChildNode] != -1) && (viableDirection(targetAtom, currentAtom, rightChildAtom, treeLevel, maxDistance) == true))
			rangeSearchCpu(targetAtom, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, rightChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);
	}


}

int squaredDistance3D(atomCoords targetAtom, atomCoords currentAtom)
{
	int temp = pow((targetAtom.xyz[0] - currentAtom.xyz[0]), 2) + pow((targetAtom.xyz[1] - currentAtom.xyz[1]), 2) + pow((targetAtom.xyz[2] - currentAtom.xyz[2]), 2);
	return temp;
}

bool viableDirection(atomCoords targetAtom, atomCoords currentAtom, atomCoords childAtom, int treeLevel, int requiredDistance)
{
	int dim = treeLevel % 3;

	if (targetAtom.xyz[dim] > currentAtom.xyz[dim])
		if ((currentAtom.xyz[dim] > childAtom.xyz[dim]) && (abs(currentAtom.xyz[dim] - targetAtom.xyz[dim])>requiredDistance))
			return false;
		else
			return true;
	else
		if ((currentAtom.xyz[dim] < childAtom.xyz[dim]) && (abs(currentAtom.xyz[dim] - targetAtom.xyz[dim])>requiredDistance))
			return false;
		else
			return true;

	return false; //this should never be hit
}
