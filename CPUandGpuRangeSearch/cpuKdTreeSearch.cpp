#include "cpuKdTreeSearch.h"

void cpuKdTreeRangeSearch(rangeSearchSettings& settings, ProteinDataHandler heldProteinSets, AtomToNumHashTable atomReferenceTable)
{

	int outputType = settings.resultsPrintFormat;
	outputHandler filePrinter;
	std::string printType;

	if (outputType == 3)
	{
		printType = "_Summary";
	}
	else if (outputType == 4)
	{
		printType = "_Detailed";
	}

	std::cout << "PERFORMING KD TREE RANGE SEARCH" << std::endl;


	int soughtAtomANumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeOne);
	int soughtAtomBNumber = atomReferenceTable.retrieveHashValue(settings.AtomTypeTwo);
	int maxDistanceSquared = settings.requiredProximity*settings.requiredProximity;


	rangeSearchArrays rangeSearch;

	int allocatedSearchSize = 16384 * sizeof(int) * 4;

	rangeSearch.atomAPositionList = (int*)malloc(sizeof(int) * 16384);
	rangeSearch.atomACount = (int*)malloc(sizeof(int) * 1);
	rangeSearch.completionFlag = (int*)malloc(sizeof(int) * 1);
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
		if (outputType == 3 || outputType == 4)
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
			clock_t n, m;
			n = clock();



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

						rangeSearchCpu(atomACoords, soughtAtomBNumber, settings.requiredProximity, resultsVector, heldProteinSets.ProteinDataHolder[i].namesSets, heldProteinSets.ProteinDataHolder[i].xCoordsSets, heldProteinSets.ProteinDataHolder[i].yCoordsSets, heldProteinSets.ProteinDataHolder[i].zCoordsSets, heldProteinSets.ProteinDataHolder[i].kdTrees, heldProteinSets.ProteinDataHolder[i].KdTreeSize, 0, 0, runCount, heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry, heldProteinSets.ProteinDataHolder[i].KdTreeSize*currentEntry);


						for (int j = 0; j < resultsVector.size(); j++) //All atomB will be related to the atomA being searched for in the greater loop
						{
							resultsSize = resultsSize + 1;


							if (outputType == 2)
							{
								std::cout << "Atom A: " << heldProteinSets.ProteinDataHolder[i].namesSets[rangeSearch.atomAPositionList[p]] << "\t - Pos : " << (rangeSearch.atomAPositionList[p] - heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry) << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000) << std::endl;
								std::cout << "Atom B: " << heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]] << "\t - Pos : " << resultsVector[j] << "\t X: " << ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000) << "\t Y: " << ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000) << "\t Z: " << ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000) << std::endl << std::endl;

							}
							if (outputType == 4)
							{
								filePrinter.printLineToOutputFile("Atom A: ", heldProteinSets.ProteinDataHolder[i].namesSets[rangeSearch.atomAPositionList[p]], "\t - Pos : ", (rangeSearch.atomAPositionList[p] - heldProteinSets.ProteinDataHolder[i].MaxEntrySize*currentEntry), "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[rangeSearch.atomAPositionList[p]])) / 1000));
								filePrinter.printLineToOutputFile("Atom B: ", heldProteinSets.ProteinDataHolder[i].namesSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]], "\t - Pos : ", resultsVector[j], "\t X: ", ((double(heldProteinSets.ProteinDataHolder[i].xCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000), "\t Y: ", ((double(heldProteinSets.ProteinDataHolder[i].yCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000), "\t Z: ", ((double(heldProteinSets.ProteinDataHolder[i].zCoordsSets[heldProteinSets.ProteinDataHolder[i].MaxEntrySize * currentEntry + resultsVector[j]])) / 1000));
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
						std::cout << "Number of matches in file " << numberOfProcessedFiles << " in set " << i << "  is: 0" << std::endl;
					else if (outputType == 3 || outputType == 4)
						filePrinter.printLineToOutputFile("Number of matches in file ", numberOfProcessedFiles, " in set ", i, "  is: 0");
					numberOfProcessedFiles++;
				}
			}

			m = clock();
			if (settings.debugLevel>0)
				print_elapsed(n, m, "time to process file subset: ");
		}
		filePrinter.closeOpenFile();
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
		int posOfRightChild = kdArray[kdTreeStartPositionOffset + rightChildNode];
		setAtomCoords(rightChildAtom, xValueArray[valueArraysStartPositionOffset + posOfRightChild], yValueArray[valueArraysStartPositionOffset + posOfRightChild], zValueArray[valueArraysStartPositionOffset + posOfRightChild]);

		int posOfLeftChild = kdArray[kdTreeStartPositionOffset + leftChildNode];
		setAtomCoords(leftChildAtom, xValueArray[valueArraysStartPositionOffset + posOfLeftChild], yValueArray[valueArraysStartPositionOffset + posOfLeftChild], zValueArray[valueArraysStartPositionOffset + posOfLeftChild]);


		if ((leftChildNode < kdArraySize) && (posOfLeftChild != -1) && (viableDirection(targetAtom, currentAtom, leftChildAtom, treeLevel, maxDistance) == true))
			rangeSearchCpu(targetAtom, elementTwo, maxDistance, resultsVector, names, xValueArray, yValueArray, zValueArray, kdArray, kdArraySize, leftChildNode, (treeLevel + 1), runcount, valueArraysStartPositionOffset, kdTreeStartPositionOffset);

		if ((rightChildNode < kdArraySize) && (posOfRightChild != -1) && (viableDirection(targetAtom, currentAtom, rightChildAtom, treeLevel, maxDistance) == true))
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
